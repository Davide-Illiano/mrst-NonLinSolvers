classdef RichardsTransportEquationFixedPointSchemes < ReservoirModel
 %% Function used to solve the coupled probelm Richards and Transport equation. 
 %  If Picard set true we use Modified Picard Method
 %  RichardsTransportEquationModelPicardScheme(G, rock, fluid, 'Picard', 1);
 %
 %  To apply the L scheme specify the 2 constant 'L_p' and 'L_c' 
 %  RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'L_p', x, 'L_c', y)
    
    properties
        % Polymer present
        polymer
        
        forces
        Picard
        Newton
        LScheme
        L_p
        L_c
        Mixed
        n
        Anderson
    end

    
    methods
        function model = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, varargin) % no varagin
            model = model@ReservoirModel(G, rock, fluid);
            model.water = true;
            model.oil = false;
            model.gas = false;
          % This is the model parameters for oil/water/polymer
            model.polymer = true;
            
            model = merge_options(model,varargin{:}); 
        end
        
        function forces = getValidDrivingForces(model)
        forces = getValidDrivingForces@ReservoirModel(model);
        %% Analytical exteranl forces
        forces.analyticalForce1 = [];
        forces.analyticalForce2 = [];
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            if model.Newton
             [problem, state] = equationsRichardsTransport(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});   
            end
            
            if model.LScheme
                [problem, state] = equationsRichardsTransportLScheme(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});  
            end
        end
    
        
            function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
        % Perform a step that ideally brings the state closer to convergence
        %
        % SYNOPSIS:
        %   [state, report] = model.stepFunction(state, state0, dt, ...
        %                                        forces, ls, nls, it)
        %
        % DESCRIPTION:
        %   Perform a single nonlinear step and report the progress towards
        %   convergence. The exact semantics of a nonlinear step varies
        %   from model to model, but the default behavior is to linearize
        %   the equations and solve a single step of the Newton-Rapshon
        %   algorithm for general nonlinear equations.

        %
        % PARAMETERS:
        %   model         - Class instance
        %   state         - Current state to be solved for time t + dt.
        %   state0        - The converged state at time t.
        %   dt            - The scalar time-step.
        %   drivingForces - Forces struct. See `getDrivingForces`.
        %   linsolver     - `LinearSolverAD` instance used to solve the
        %                   linear systems that may appear from
        %                   linearization.
        %   nonlinsolver  - `NonLinearSolverAD` controlling the solution
        %                   process.
        %   iteration     - The current nonlinear iterations number. Some
        %                   models implement special logic depending on the
        %                   iteration (typically doing setup during the
        %                   first iteration only).
        %
        % OPTIONAL PARAMETERS:
        %   varargin - Any additional arguments are passed onto
        %              `getEquations` without modification or validation.
        %
        % RETURNS:
        %   state  - Updated state `struct` that hopefully is closer to
        %            convergence in some sense.
        %   report - A report produced by `makeStepReport` which indicates
        %            the convergence status, residual values and other
        %            useful information from the application of the
        %            `stepFunction` as well as any dispatched calls to
        %            other functions.
        %
        % SEE ALSO:
        %   `NonLinearSolverAD`, `LinearSolverAD`, `simulateScheduleAD`
        %
        state_start = state;
        
        onlyCheckConvergence = iteration > nonlinsolver.maxIterations;
        timer = tic();
        [problem, state] = model.getEquations(state0, state, dt, drivingForces, ...
                                   'ResOnly', onlyCheckConvergence, ...
                                   'iteration', iteration, ...
                                   varargin{:});
                               
        problem.iterationNo = iteration;
        problem.drivingForces = drivingForces;
        t_assembly = toc(timer);
        
        %% We add now the analytical external forces directly to the equations
        problem.equations{1} = problem.equations{1}-dt.*drivingForces.analyticalForce1;
        problem.equations{2} = problem.equations{2}-dt.*drivingForces.analyticalForce2;
        
        [convergence, values, resnames] = model.checkConvergence(problem);
        
        % Minimum number of iterations can be prescribed, i.e., we
        % always want at least one set of updates regardless of
        % convergence criterion.
        doneMinIts = iteration > nonlinsolver.minIterations;

        % Defaults
        failureMsg = '';
        failure = false;
        [linearReport, updateReport, stabilizeReport] = deal(struct());
        if (~(all(convergence) && doneMinIts) && ~onlyCheckConvergence)
            % Get increments for Newton solver
            [dx, ~, linearReport] = linsolver.solveLinearProblem(problem, model);
            if any(cellfun(@(d) ~all(isfinite(d)), dx))
                failure = true;
                failureMsg = 'Linear solver produced non-finite values.';
            end
            % Let the nonlinear solver decide what to do with the
            % increments to get the best convergence
            [dx, stabilizeReport] = nonlinsolver.stabilizeNewtonIncrements(model, problem, dx);

            if (nonlinsolver.useLinesearch && nonlinsolver.convergenceIssues) || ...
                nonlinsolver.alwaysUseLinesearch
                [state, updateReport, stabilizeReport.linesearch] = nonlinsolver.applyLinesearch(model, state0, state, problem, dx, drivingForces, varargin{:});
            else
                % Finally update the state. The physical model knows which
                % properties are actually physically reasonable.
                if model.Anderson == 1 
                    linearReport.state_prev = state;
                end
                [state, updateReport] = model.updateState(state, problem, dx, drivingForces);
            end
        end
        isConverged = (all(convergence) && doneMinIts) || model.stepFunctionIsLinear;
        
        % If step function is linear, we need to call a residual-only
        % equation assembly to ensure that indirect/derived quantities are
        % set with the updated values (fluxes, mobilities and so on).
        if model.stepFunctionIsLinear
            [~, state] = model.getEquations(state0, state, dt, drivingForces, ...
                                   'ResOnly', true, ...
                                   'iteration', iteration+1, ...
                                   varargin{:});
        end
        if model.verbose
            printConvergenceReport(resnames, values, convergence, iteration);
        end
        %% Necessary for richards+transport
        %norm_p = norm(state_start.pressure-state.pressure);
        %norm_c = norm(state_start.c-state.c);
        
        report = model.makeStepReport(...
                        'LinearSolver', linearReport, ...
                        'UpdateState',  updateReport, ...
                        'Failure',      failure, ...
                        'FailureMsg',   failureMsg, ...
                        'Converged',    isConverged, ...
                        'AssemblyTime', t_assembly, ...
                        'Residuals',    values, ...
                        'StabilizeReport', stabilizeReport,...
                        'ResidualsConverged', convergence);
                    
    end

        
         function state = validateState(model, state)
            state = validateState@ReservoirModel(model, state);
            % Polymer must be present
            model.checkProperty(state, 'concentration', model.G.cells.num, 1);
        end

        function [state, report] = updateState(model, state, problem, ...
                dx, drivingForces)
            [state, report] = updateState@ReservoirModel(model, ...
               state, problem,  dx, drivingForces);
           if model.polymer
                c = model.getProp(state, 'c');
            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ReservoirModel(model, state0, state, dt, drivingForces);
            if model.polymer
                c     = model.getProp(state, 'c');
            end
        end       
       
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            switch(lower(name))
                case {'c', 'concentration'} %, 'polymer'
                    fn = 'c';
                    index = 1;
                case {'cmax'}
                    index = 1;
                    fn = 'cmax';
               % case 'qwc'
               %     index = 1;
               %     fn = 'qWc';
                otherwise
                    
                    [fn, index] = getVariableField@ReservoirModel(model, name);
            end
        end
        
        function names = getComponentNames(model)
            names = getComponentNames@ReservoirModel(model);
            if model.polymer
                names{end+1} = 'concentration';
            end
        end
       
       function [names, types] = getExtraWellEquationNames(model)
            [names, types] = getExtraWellEquationNames@ReservoirModel(model);
            if model.polymer
                names{end+1} = 'cWells';
                types{end+1} = 'perf';
            end
       end
       
       
       function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@ReservoirModel(model);
            if model.polymer
                names{end+1} = 'qWc';
            end
       end
      
      
       function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ReservoirModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.polymer
                assert(model.water, 'Polymer injection requires a water phase.');
                f = model.fluid;

                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition

                if well.isInjector
                    concWell = model.getProp(well.W, 'c');
                    cqP = concWell.*cqWs;
                else
                    pix = strcmpi(model.getComponentNames(), 'c');
                    concWell = packed.components{pix};

                    a = f.muWMult(f.cmax).^(1-f.mixPar);
                    cbarw     = concWell/f.cmax;

                    % the term (a + (1 - a).*cbarw) account for the
                    % todd-longstaff mixing factor, which model the fact that for
                    % not-fully mixed polymer solution the polymer does not
                    % travel at the same velocity as water. See the governing
                    % equation for polymer (e.g. equationsOilWaterPolymer.m)
                    cqP = concWell.*cqWs./(a + (1-a).*cbarw);
                end

                qwpoly = packed.extravars{strcmpi(packed.extravars_names, 'qwc')};

                compEqs{end+1} = qwpoly - sum(concWell.*cqWs);
                compSrc{end+1} = cqP;
                eqNames{end+1} = 'cWells';
            end
        end
        %}
        %
            function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % INPUT:
        %
        % model  - (Base class, automatic)
        %
        % cname  - Name of the component. Must be a property known to the
        %          model itself through getProp/getVariableField.
        %
        % eq     - Equation where the source terms are to be added. Should
        %          be one value per cell in the simulation grid (model.G)
        %          so that the src.sourceCells is meaningful.
        %
        % component - Cell-wise values of the component in question. Used
        %          for outflow source terms only.
        %
        % src    - Source struct containing fields for fluxes etc. Should
        %          be constructed from force and the current reservoir
        %          state by computeSourcesAndBoundaryConditionsAD.
        %
        % force  - Force struct used to produce src. Should contain the
        %          field defining the component in question, so that the
        %          inflow of the component through the boundary condition
        %          or source terms can accurately by estimated.
        if isempty(force)
            return
        end
        c = model.getProp(force, cname);
        cells = src.sourceCells;
        switch lower(cname)
            case {'concentration'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
            otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
        end
        eq(cells) = eq(cells) - qC;
        src.components{end+1} = qC;
    end
     %} 
     
        function p = getProp(model, state, name)
        % Get a single property from the nonlinear state
        %
        % SYNOPSIS:
        %   p = model.getProp(state, 'pressure');
        %
        % PARAMETERS:
        %   model - Class instance.
        %   state - `struct` holding the state of the nonlinear problem.
        %   name  - A property name supported by the model's
        %           `getVariableField` mapping.
        %
        % RETURNS:
        %   p     - Property taken from the state.
        %
        % SEE ALSO:
        %   `getProps`
        
        [fn, index] = model.getVariableField(name);
        p = state.(fn)(:, index);
        end
    
    end
  
 

    
end