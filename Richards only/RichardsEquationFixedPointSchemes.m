classdef RichardsEquationFixedPointSchemes < ReservoirModel
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
        Newton
        LScheme
        L_p
        L_c
        Mixed
        n
        analyticalForce1
        Anderson
    end

    
    methods
        function model = RichardsEquationFixedPointSchemes(G, rock, fluid, varargin) % no varagin
            model = model@ReservoirModel(G, rock, fluid);
            model.water = true;
            model.oil = false;
            model.gas = false;
          % This is the model parameters for oil/water/polymer (used in case of transport equation)
            model.polymer = false;
            
            model = merge_options(model,varargin{:}); 
        end
        
        function forces = getValidDrivingForces(model)
        forces = getValidDrivingForces@ReservoirModel(model);
        
        %% Analytical exteranl forces
        forces.analyticalForce1 = [];
        end
        
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

          
            if model.Newton
             [problem, state] = equationsRichards(state0, state, model,...
                            dt, ...
                            drivingForces,...
                            varargin{:});   
            end   
            if model.LScheme
                [problem, state] = equationsRichardsLScheme(state0, state, model,...
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

        end
  
end