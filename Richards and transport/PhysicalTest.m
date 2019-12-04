mrstModule add ad-blackoil ad-core ad-props mrst-gui %blackoil-sequential diagnostics
%close all
clear all
clc
mrstVerbose off
gravity on


% Create cartesian grid
for J = [10]
X = J; % mesh size
G = cartGrid([X X], [100 100]);  %2D problem easily extendable to 3D
G = computeGeometry(G);

rho = 1; 
mu = 0.01;                  % [g/cm.s] water viscosity
g = 980.66; 
K_sat = 33.1920/hour;      % [cm/s] saturated hydraulic conductivity
k = (K_sat*mu)/(rho*g);

x = G.cells.centroids(:,1);
y = G.cells.centroids(:,2);

[ii, jj] = gridLogicalIndices(G);
lower = jj < X*1/4;  % vadose zone
higher = jj > X*1/4;
middle = jj == X*1/4;

K = k*ones(G.cells.num, 1);
%k(lower) = y(lower).*0 + k .* 1e3;



%poro = gaussianField(G.cartDims, [.2 .6],[11 3], 3.5);
%K = poro.^3.*(1e-6)./(.81*72*(1-poro).^2);
%rock = makeRock(G, K(:),poro(:));

rock = struct('perm', K(:), ...
              'poro', .5*ones(G.cells.num, 1)); 

%max(max(K)/min(min(K)))
%{
% show poro
figure;
plotCellData(G,rock.poro,'EdgeColor', 'none');
colorbar; axis equal tight;
title('Rock porosity')
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',16,...
'FontName','Times')

%show perm
figure;
plotCellData(G,rock.perm);
colorbar; axis equal tight;
title('Rock permeability');
set(gca,...
'Units','normalized',...
'FontUnits','points',...
'FontWeight','normal',...
'FontSize',16,...
'FontName','Times')
%}

% Default water fluid with unit values
fluid = initSimpleADIFluid('phases', 'W');

% Define the water content theta and Kmult as function of the pressure 
fluid.theta = @(p, c) getThetaCoupled(p, c, 0.102, 0.368, 0.0335, 2, .44, .0046);  % (p, c, theta_R, theta_S, alpha,  n, a,b)
fluid.Kmult = @(p, theta) getConductivity(p, theta, 1, 2);   %.8 is K_s probably unrealistic number, required in def of K

%% Fixed point solver, choose linearization scheme: Newton or LScheme
disp('We start by investigating the unsaturated domain.')
choice = menu('Which linearization scheme do you want to use? Please choose one','Newton','L-Scheme');

if choice==1 
    
    disp('You choose to use the Newton method.')
    modelFP = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'Newton', 1 );

elseif choice==2

    disp('You choose to use the L-Scheme method.')
    modelFP = RichardsTransportEquationFixedPointSchemes(G, rock, fluid, 'LScheme', 1 );
    modelFP.L_p = .005;
    modelFP.L_c = .015;
    
else
    disp('Not valid choice')
    return;
end

state0 = initResSol(G, -y + 10);      % initial pressure p and c
state0.c = 1.*ones(G.cells.num,1);

%% Time domain and time step
time = 72*hour;
n = 10;
dT = time/n;

%% Boundary conditions
bc = [];
bc = fluxside(bc, G, 'ymax', 1.*meter^3/hour(), 'sat', 1);
bc.c = 10*ones(size(bc.sat,1), 1);

bc = fluxside(bc, G, 'ymin', -1.*meter^3/hour(), 'sat', 1);
bc.c = 1*ones(size(bc.sat,1), 1);

schedule = simpleSchedule(repmat(dT,1,n),'bc', bc);

modelFP.nonlinearTolerance = 1e-5; 
[modelFP.forces] = getValidDrivingForces(modelFP);


%% Introduce external forces
for i=1:n
   analyticalForce1 = 0;   %we want zero external forces    
   analyticalForce2 = 0;
   
   modelFP.forces(i).analyticalForce1 = analyticalForce1; % original: .006.*cos(4/3*pi.*y).*sin(x);
   modelFP.forces(i).analyticalForce2 = analyticalForce2;
end

nls = NonLinearSolver('maxIterations', 10000,'enforceResidualDecrease', true);
%
tic
[~,statesFP, reportFP] = FixedPointSimulateScheduleAD(state0, modelFP, schedule, 'nonlinearsolver', nls);
toc

disp('Total number iterations:')
sum(reportFP.Iterations)
%
%
%% Plot the result
mrstModule add mrst-gui
close all
plotToolbar(G, statesFP)
title('Newton')
colorbar
%}

modelFP.Anderson = 1;
if modelFP.Anderson == 1
    nls = NonLinearSolverAnderson('maxIterations', 2000,'enforceResidualDecrease', true);
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

end