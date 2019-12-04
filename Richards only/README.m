%% Flow in variably saturated media
% Here we investigate the flow of water in unsaturated/saturated media,
% described by the Richars equation.
% Such equation degenerates whenever we enter in teh fully saturated part
% of the doamin. We will observe that in these cases the newton method will
% fail to converge.

% Files included:
% 1) Example_variablysaturated_Domain.m : main file where the problem is defined.
% We are asked to choose which linearization scheme (Newton/L-scheme) we want
% to implement and if, after having solved the problem on the unsaturated
% domain, we want to exted the study to the variably saturated domain.
% 
% 2) equationsRichards.m: the richards equation is here defined and the
% Newton method is used.
% 
% 3) equationsRichardsLScheme: we present the linearized Richards equation
% obtained thanks to the LScheme
% 
% 4) FixedPointSimulateScheduleAD.m: we modified the origina
% SimulateScheduleAD.m to include analythical external forces used to
% define our problem.
% 
% 5) getTheta.m and getConductivity.m: two files to obtain the water
% content theta, expressed as a function of the unknown pressure, and the
% conductivity K, espressed as function of both pressure and theta. The
% equations are based on the van Genuchten formulation.
% 
% 6) PlotsIterationsResiduals.m: file used to plot the iterations and the
% residual obtained into Example_variablysaturated_Domain.m.
% 
% 7) NonLinearSolverAnderson.m: new file obtained modifying the original NonLinearSolver.m 
% Adding a couple of lines we could introduce the loop used to obtain the
% Anderson acceleration.
