classdef LinearSolverAD < handle
    % Base class for linear solvers in the AD framework
    %
    % SYNOPSIS:
    %   solver = LinearSolverAD()
    %
    % DESCRIPTION:
    %   Base class for linear solvers. Implement methods for solving linearized
    %   problems and adjoints. Also supports setup/cleanup functions
    %   before/after use for initialize once-type usage.
    %
    % REQUIRED PARAMETERS:
    %   None
    %
    % OPTIONAL PARAMETERS:
    %   'property' - Set property to value.
    %
    % NOTE:
    %   This class is intended as superclass. It cannot actually solve
    %   problems.
    %
    % SEE ALSO:
    %   `BackslashSolverAD`, `CPRSolverAD`, `LinearizedProblem`
    
    properties
        tolerance % Linear solver tolerance
        maxIterations % Max number of iterations used
        extraReport % Enable this to produce additional report output
        % May use a lot of memory for large problems
        verbose % Verbose output enabler
        replaceNaN % Boolean indicating if the solver should replace NaN in the results
        replaceInf % Boolean indicating if the solver should replace Inf in the results
        replacementNaN % If replaceNaN is enabled, this is the value that will be inserted
        replacementInf % If replaceInf is enabled, this is the value that will be inserted
        reduceToCell % Reduce to per-cell system before solving
        applyLeftDiagonalScaling % Apply left diagonal scaling before solving
        applyRightDiagonalScaling % Apply right diagonal scaling before solving
        keepNumber % If set, linear solver will reduce the system to the first keepNumber entries
    end
    
    methods
        function solver = LinearSolverAD(varargin)
            solver.tolerance       = 1e-8;
            solver.maxIterations   = 100000; %25
            solver.extraReport     = false;
            solver.verbose         = mrstVerbose();
            solver.replaceNaN      = false;
            solver.replaceInf      = false;
            solver.replacementNaN  = 0;
            solver.replacementInf  = 0;
            solver.reduceToCell    = false;
            solver.applyLeftDiagonalScaling = false;
            solver.applyRightDiagonalScaling = false;
            
            solver = merge_options(solver, varargin{:});
            
            assert(solver.maxIterations >= 0);
            assert(solver.tolerance >= 0);
        end
        
        function [result, report] = solveLinearSystem(solver, A, b) %#ok
            % Solve the linear system to a given tolerance
            error('Superclass not meant for direct use')
        end
        
        function [grad, result, report] = solveAdjointProblem(solver, problemPrev,...
                problemCurr, adjVec, objective, model) %#ok
            % Solve an adjoint problem.
            problemCurr = problemCurr.assembleSystem();
            
            % Maybe formalize the control variables a bit in the future
            % sometime...
            if iscell(objective)
                objective = objective{:};
            end
            objective = combineEquations(objective);
            assert(isa(objective, 'ADI'), 'Objective function was not of type ADI.');
            rhs = -(objective.jac{1})';
            if ~isempty(adjVec)
                problemPrev = problemPrev.assembleSystem();
                rhs = rhs - problemPrev.A'*adjVec;
            end
            A = problemCurr.A';
            [result, report] = solver.solveLinearSystem(A, rhs);
            grad = solver.storeIncrements(problemCurr, result);
        end
        
        function [dx, result, report] = solveLinearProblem(solver, problem, model)
            % Solve a linearized problem
            if solver.reduceToCell
                % Eliminate non-cell variables (well equations etc)
                keep = problem.indexOfType('cell');
                [problem, eliminated] = solver.reduceToVariable(problem, keep);
            end
            problem = problem.assembleSystem();
            assert(all(isfinite(problem.b)), 'Linear system rhs must have finite entries.');
            
            timer = tic();
            % Get linearized system
            [A, b] = problem.getLinearSystem();
            % Disp condition number of system
            %{
            disp('Condition number:')
            condest(A)
            pause
            %}
            % Apply scaling
            [A, b, scaling] = solver.applyScaling(A, b);
            % Reduce system (if requested)
            [A, b, lsys] = solver.reduceLinearSystem(A, b);
            % Solve the system
            [result, report] = solver.solveLinearSystem(A, b); 
            %norm(result)
            %pause
            % Recover eliminated variables on linear level
            result = solver.recoverLinearSystem(result, lsys);
            % Undo scaling
            result = solver.undoScaling(result, scaling);
            
            [result, report] = problem.processResultAfterSolve(result, report);
            report.SolverTime = toc(timer);
            if solver.replaceNaN
                result(isnan(result)) = solver.replacementNaN;
            end
            if solver.replaceInf
                result(isinf(result)) = solver.replacementInf;
            end
            dx = solver.storeIncrements(problem, result);
            if solver.reduceToCell
                dx = solver.recoverResult(dx, eliminated, keep);
            end
            %% Store residuals for anderson Acceleration
            report.VectorResiduals = dx;
        end
        
        function dx = storeIncrements(solver, problem, result) %#ok
            % Extract the results from a vector into a cell array with one
            % entry per primary variable in the linearized problem.
            
            % Find first index corresponding to ADI equation
            ix = find(cellfun(@(x) isa(x, 'ADI'), problem.equations), 1);
            if isempty(ix)
                % No ADI equations, we return empty increments
                assert(isempty(result), ...
                    'No ADI equations returned. Unable to map increments to variables.');
                dx = {};
                return
            end
            % Calculate positions in newton increment
            numVars = problem.equations{ix}.getNumVars();
            cumVars = cumsum(numVars);
            ii = [[1;cumVars(1:end-1)+1], cumVars];
            
            eqn = size(ii,1);
            dx = cell(eqn,1);
            for i = 1:eqn
                dx{i} = result(ii(i,1):ii(i,2), :);
            end
        end
        
        function [A, b, sys] = reduceLinearSystem(solver, A, b)
            % Perform Schur complement reduction of linear system
            
            sys = struct('B', [], 'C', [], 'D', [], 'f', [], 'h', [], 'E_L', [], 'E_U', []);
            if isempty(solver.keepNumber) || solver.keepNumber >= size(b, 1)
                return
            end
            if isempty(A)
                return
            end
            nk = solver.keepNumber;
            start = 1:nk;
            [ix, jx, vx] = find(A);
            if any(~isfinite(vx))
                warning('Non-finite values in matrix before Schur-complement reduction.');
            end
            if any(~isfinite(b))
                warning('Non-finite values in right-hand side before Schur-complement reduction.');
            end
            n = size(A, 2);
            keep = false(n, 1);
            keep(start) = true;
            
            
            keepRow = keep(ix);
            keepCol = keep(jx);
            kb = keepRow & keepCol;
            sys.B = sparse(ix(kb), jx(kb), vx(kb), nk, nk);
            
            kc = keepRow & ~keepCol;
            sys.C = sparse(ix(kc), jx(kc) - nk, vx(kc), nk, n - nk);
            
            kd = ~keepRow & keepCol;
            sys.D = sparse(ix(kd) - nk, jx(kd), vx(kd), n - nk, nk);
            
            ke = ~keepRow & ~keepCol;
            sys.E = sparse(ix(ke) - nk, jx(ke) - nk, vx(ke), n - nk, n - nk);
            sys.f = b(keep);
            sys.h = b(~keep);
            
            [sys.E_L, sys.E_U] = lu(sys.E);
            A = sys.B - sys.C*(sys.E_U\(sys.E_L\sys.D));
            b = sys.f - sys.C*(sys.E_U\(sys.E_L\sys.h));
        end
        
        function x = recoverLinearSystem(solver, x, sys)
            % Recover eliminated variables
            if ~isempty(sys.E_U)
                s = sys.E_U\(sys.E_L\(sys.h - sys.D*x));
                x = [x; s];
            end
        end
        
        function [A, b, scaling] = applyScaling(solver, A, b)
            % Apply left or right diagonal scaling
            scaling = struct();
            applyLeft = solver.applyLeftDiagonalScaling;
            applyRight = solver.applyRightDiagonalScaling;
            if ~applyLeft && ~applyRight
                return
            end
            M = solver.getDiagonalInverse(A);
            if solver.applyLeftDiagonalScaling
                assert(~applyRight, 'Cannot both apply left and right diagonal scaling');
                A = M*A;
                b = M*b;
            else
                A = A*M;
            end
            scaling.M = M;
        end
        
        function x = undoScaling(solver, x, scaling)
            % Undo effects of scaling applied to linear system
            if solver.applyRightDiagonalScaling
                x = scaling.M*x;
            end
        end
        
        function x = preconditionerInverse(solver, M, x)
            % Apply a preconditioner. Either a handle or a matrix.
            if isempty(M)
                return
            end
            
            if isa(M, 'function_handle')
                % We got a function handle for the inverse
                x = M(x);
            else
                % We got a matrix
                x = M\x;
            end
        end
        
        function solver = setupSolver(solver, A, b, varargin) %#ok
            % Run setup on a solver for a given system
            
            % Dummy function run before a set of linear problems with
            % different right hand sides
        end
        
        function  solver = cleanupSolver(solver, A, b, varargin) %#ok
            % Clean up solver after use (if needed)
            
            % Dummy function run after a set of linear problems with
            % different right hand sides
            
            % For overloading when for example calling multigrid solvers as
            % a preconditioner
        end
        
        function [dx, result, report] = solveCellReducedLinearProblem(solver, problem, model)
            % Reduce a problem to cell-variables, solve and then recover
            % the eliminated variables afterwards.
            
            % Eliminate non-cell variables (well equations etc)
            keep = problem.indexOfType('cell');
            
            [problem, eliminated] = solver.reduceToVariable(problem, keep);
            
            % Solve a linearized problem
            problem = problem.assembleSystem();
            
            timer = tic();
            [result, report] = solver.solveLinearSystem(problem.A, problem.b);
            report.SolverTime = toc(timer);
            
            dxCell = solver.storeIncrements(problem, result);
            
            dx = solver.recoverResult(dxCell, eliminated, keep);
        end
        
        function [problem, eliminated] = reduceToVariable(solver, problem, keep)
            remove = find(~keep);
            
            problem = problem.clearSystem();
            
            eliminated = cell(numel(remove), 1);
            elimNames = problem.equationNames(remove);
            
            for i = 1:numel(remove)
                [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
            end
        end
        
        function dx = recoverResult(solver, dxElim, eliminatedEqs, keep)
            kept = find(keep);
            left = find(~keep);
            keptEqNo = numel(kept);
            
            % Find number of variables
            nP = numel(keep);
            
            % Set up storage for all variables, including those we
            % eliminated previously
            dx = cell(nP, 1);
            
            % Recover non-cell variables
            recovered = false(nP, 1);
            recovered(kept) = true;
            
            % Put the recovered variables into place
            dx(recovered) = dxElim;
            
            for i = numel(eliminatedEqs):-1:1
                pos = left(i);
                dVal = recoverVars(eliminatedEqs{i}, keptEqNo + 1, dx(recovered));
                dx{pos} = dVal;
                recovered(pos) = true;
            end
        end
        
        function M = getDiagonalInverse(solver, A)
            % Reciprocal of diagonal matrix
            sz = size(A);
            assert(sz(1) == sz(2), 'Matrix must be square!');
            n = sz(1);
            d = 1./diag(A);
            d(~isfinite(d)) = 1;
            I = (1:n)';
            M = sparse(I, I, d, n, n);
        end
    end
end

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
