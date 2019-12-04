mrstModule add ad-blackoil ad-core ad-props mrst-gui %blackoil-sequential diagnostics
close all
clear all
clc
mrstVerbose on
gravity on


% Create cartesian grid
X = 10; % mesh size
G = cartGrid([X X], [1 1]*meter);  %2D problem easily extendable to 3D
G = computeGeometry(G);

% Homogenous rock properties
rock = struct('perm', 1000*darcy*ones(G.cells.num, 1), ...
              'poro', .5*ones(G.cells.num, 1)); 

% Default water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');

% Define the water content theta and Kmult as function of the pressure 
fluid.theta = @(p, c) getThetaCoupled(p, c, 0.026, 0.42, 0.551, 2.9, .44, .0046);  % (p, c, theta_R, theta_S, alpha,  n, a,b)
fluid.Kmult = @(p, theta) getConductivity(p, theta, .12, 2.9);   %.8 is K_s probably unrealistic number, required in def of K

%% Fixed point solver, choose linearization scheme: Newton or LScheme
disp('We start by investigating the unsaturated domain.')
choice = menu('Which linearization scheme do you want to use? Please choose one','Newton','L-Scheme');

if choice==1 
    
    disp('You choose to use the Newton method.')
    modelFP = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'Newton', 1 );

elseif choice==2

    disp('You choose to use the Newton method.')
    modelFP = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'LScheme', 1 );
    modelFP.L_p = .15;
    modelFP.L_c = .15;
    
else
    disp('Not valid choice')
    return;
end

%% Compute external forces, f(x,y) and iitial pressure
x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

[ii, jj] = gridLogicalIndices(G);

lower = jj <= X*1/4;  
higher = jj > X*1/4;  
p = -2*ones(G.cells.num, 1);
p(lower) = -y(lower) - 1/4; % produce unsaturated domain, the pressure remains negative the domain

p_var = -2*ones(G.cells.num, 1);
p_var(lower) = -y(lower) + 1/4; % produce variably saturated domain,
                                % the pressure becomes positive in the
                                % lower part of the domain

state0 = initResSol(G, p);      % initial pressure p and c(unsaturated problem)
state0.c = 0*ones(G.cells.num,1);
state0_var = initResSol(G, p_var); % initial pressure p and c(variably saturated problem)
state0_var.c = 0*ones(G.cells.num,1);

%% Time domain and time step
time = 1;
n = 20;
dT = time/n;

%% Boundary conditions
bc = [];
bc = pside(bc, G, 'ymax', -2, 'sat', 1);%-3
bc.c = 4*ones(size(bc.sat,1), 1);
schedule = simpleSchedule(repmat(dT,1,n),'bc', bc);

modelFP.nonlinearTolerance = 1e-6; 
[modelFP.forces] = getValidDrivingForces(modelFP);


%% Introduce external forces
for i=1:n
   analyticalForce1 = .006.*cos(-4/3*pi.*y).*sin(x);   % we want external forces in vadose zone
   analyticalForce1(lower) = 0;                        % lower = jj > X*1/4;
   
   analyticalForce2 = .006.*cos(-4/3*pi.*y).*sin(x);
   analyticalForce2(lower) = 0;         
   
   modelFP.forces(i).analyticalForce1 = analyticalForce1; 
   modelFP.forces(i).analyticalForce2 = analyticalForce2;
end

nls = NonLinearSolver('maxIterations', 20,'enforceResidualDecrease', true);
%
tic
[~,statesFP, reportFP] = FixedPointSimulateScheduleAD(state0, modelFP, schedule, 'nonlinearsolver', nls);
toc

disp('Total number iterations:')
sum(reportFP.Iterations)
%}
%{
%% Plot the result
%% Plot the result Newton
mrstModule add mrst-gui
close all
plotToolbar(G, states)
title('Newton')
colorbar
%}

%% Move now to second problem given by the variably saturated domain
choice = menu('Do you want to proceed modeling the variably saturated domain? Press yes or no','Yes','No');
if choice==2 || choice==0
    disp('You choose to not investigate the variably saturated domain')
   return;
else
    disp('Solving the varuably saturated domain')
    
    choice = menu('Which linearization scheme do you want to use? Please choose one','Newton','L-Scheme');
    if choice==1 
    
        disp('You choose to use the Newton method.')
        modelFPvar = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'Newton', 1 );    %'Mixed', 1, 'L_p', .2, 'L_c', .01 

    elseif choice==2

        disp('You choose to use the L-Scheme.')
        modelFPvar = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'LScheme', 1 );    %'Mixed', 1, 'L_p', .2, 'L_c', .01 
        modelFPvar.L_p = .15;
        modelFPvar.L_c = .005;
    
    end
end

modelFPvar.nonlinearTolerance = modelFP.nonlinearTolerance; 
[modelFPvar.forces] = getValidDrivingForces(modelFPvar);

for i=1:n
   analyticalForce1 = .006.*cos(-4/3*pi.*y).*sin(x);    % we want external forces in vadose zone
   analyticalForce1(lower) = 0;                         % lower = jj > X*1/4;
   
   analyticalForce2 = .006.*cos(-4/3*pi.*y).*sin(x);
   analyticalForce2(lower) = 0;      
   
   modelFPvar.forces(i).analyticalForce1 = analyticalForce1; 
   modelFPvar.forces(i).analyticalForce2 = analyticalForce2;
   
end
%
[~,statesFPvar, reportFPvar] = FixedPointSimulateScheduleAD(state0_var, modelFPvar, schedule, 'nonlinearsolver', nls);
disp('Total number iterations:')
sum(reportFPvar.Iterations)
%}
%{
disp('Variably saturated problem has been solved')
disp('Total number iterations for unsaturated problem is:')
sum(reportFPvar.Iterations)
%}

modelFP.Anderson = 1;
modelFPvar.Anderson = 1;
if modelFP.Anderson == 1 || modelFPvar.Anderson == 1
    nls = NonLinearSolverAnderson('maxIterations', 20,'enforceResidualDecrease', true);
end
%
choice = menu('Do you want to solve once more the unsaturated domain using the L-Scheme together with the Anderson acceleration? Press yes or no','Yes','No');
if choice==2 || choice==0
    disp('You choose to not investigate the problem')
   return;
else

[~,statesFP, reportFP] = FixedPointSimulateScheduleAD(state0, modelFP, schedule, 'nonlinearsolver', nls);
disp('Unsaturated problem has been solved')
disp('Total number iterations for unsaturated problem solved with L-Scheme and Anderson acceleration:')
sum(reportFP.Iterations)
end
%}

choice = menu('Do you want to solve once more the variably saturated domain using the L-Scheme together with the Anderson acceleration? Press yes or no','Yes','No');
if choice==2 || choice==0
    disp('You choose to not investigate the problem')
   return;
else
[~,statesFPvar, reportFPvar] = FixedPointSimulateScheduleAD(state0_var, modelFPvar, schedule, 'nonlinearsolver', nls);
disp('Variably saturated problem has been solved')
disp('Total number iterations for unsaturated problem solved with L-Scheme and Anderson acceleration:')
sum(reportFPvar.Iterations)
end