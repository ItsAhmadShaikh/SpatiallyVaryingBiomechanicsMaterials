%%

%open DEMO_febio_0001_cube_uniaxial.m
%open DEMO_febio_0042_inverse_FEA_cube_uniaxial.m
%(very optional) open DEMO_febio_0027_layer_spatially_varying_material.m

%%
% % % % This is a test for Neo Hookean Material Properties
clear; close all; clc;
deleteFile = 1;
sampleDim = [1.5 1.5 1.5]*1e-3; %x y z
sampleMesh = [10 10 10 1]; %x y z_big_layers z_per_big_layer
E_ini = [1.174 3.522 2.348]*1e6; % E
thickness_ini = [0 0.675 1.5]*1e-3; % thickness

% E_ini = [1.174 1.174 3.522 3.522 2.348 2.348]; % E
% thickness_ini = [0 0.3 0.3+0.005 0.3+0.005+0.75 0.3+.0005+0.75+0.005 1.5]; % thickness

BCtype= 1; %1- symmetry, 2- face fixed Use onl1 1 if ratio =/= 0
strain = 0.05;
ratio = 0.75; %between X and Y dy/dx
rho_mat = 1.12e3;

[X2,Y2]=UniaxialTestFunction(deleteFile,sampleDim,sampleMesh, E_ini, thickness_ini, BCtype,strain,ratio, rho_mat);

%%
% % % % This is a test for Neo Hookean Material Properties with half the strain
% clear;
close all; 
clc;
deleteFile = 1;
sampleDim = [15 15 1.5]*1e-3; %x y z
sampleMesh = [10 10 10 10]; %x y z_big_layers z_per_big_layer
% E_ini = [1.174 3.522 2.348]; % E
% thickness_ini = [0 0.675 1.5]; % thickness

E_ini = [1.174 1.174 3.522 3.522 2.348 2.348]*1e6; % E
thickness_ini = [0 0.3 0.3+0.005 0.3+0.005+0.75 0.3+.0005+0.75+0.005 1.5]*1e-3; % thickness

strain = 0.1;
ratio = 1; %between X and Y dy/dx
rho_mat = 1.12e3; %1.12e-9;

[X2,Y2]=UniaxialTestFunction_Half(deleteFile,sampleDim,sampleMesh, E_ini, thickness_ini,strain,ratio, rho_mat);

%%

% % % % This is a test for OGH (Gasser 2006) Material Properties

% <div class="csl-entry">Amabili, M., Balasubramanian, P., Bozzo, I., Breslavsky, I. D., &#38; Ferrari, G. (2019). 
% Layer-specific hyperelastic and viscoelastic characterization of human descending thoracic aortas. 
% <i>Journal of the Mechanical Behavior of Biomedical Materials</i>, <i>99</i>, 27–46. 
% https://doi.org/10.1016/j.jmbbm.2019.07.008</div>

% <material id="2" type="Holzapfel-Gasser-Ogden">
%   <c>7.64</c>
%   <k1>996.6</k1>
%   <k2>524.6</k2>
%   <gamma>49.98</gamma>
%   <kappa>0.226</kappa>
%   <k>1e5</k>
% </material>

clear; close all; clc;

%Donor 3 thickness, avg prop
intima = 0.29e-3;
media = 1.13e-3;
adventitia = 0.75e-3;

c = [37.5 47.3 23.3]*1e3*2;
k1 = [11.6 0.448 4.08]*1e6*2;
k2 = [107.5 35.19 87.4];
gamma = [0.809 0.859 0.728];
kappa = [0.17 0.24 0.24];
k = [435415260 435415260 435415260];

deleteFile = 1;
sampleDim = [intima+media+adventitia intima+media+adventitia intima+media+adventitia]; %x y z
sampleMesh = [10 10 10 1]; %x y z_big_layers z_per_big_layer
% E_ini = [1.174 3.522 2.348]; % E
% thickness_ini = [0 0.675 1.5]; % thickness
thickness_ini = [0 intima+media/2 intima+media+adventitia]; % thickness


BCtype= 1; %1- symmetry, 2- face fixed Use onl1 1 if ratio =/= 0
strain = 0.2;
ratio = 1; %between X and Y dy/dx
rho_mat = 1.12e-9;
[X2,Y2]=UniaxialTestFunctionOGH(deleteFile,sampleDim,sampleMesh, c, k1, k2, gamma, kappa, k, thickness_ini, BCtype,strain,ratio, rho_mat);




%%
% This sectiton is the first test of optimizaition

clear; close all; clc;
deleteFile = 0;
sampleDim = [1.5 1.5 1.5]; %x y z
sampleMesh = [10 10 10 1]; %x y z_big_layers z_per_big_layer
E_ini1 = [1.174 3.522 2.348]; % E
thickness_ini = [0 0.675 1.5]; % thickness

ini_guess(:,1) = [1.174 3.522 2.348]; % E
ini_guess(:,2) = [0 0.675 1.5]; % thickness

global t;
t = 0;


% [X,Y]=Test5_Optimization(ini_guess,deleteFile,sampleDim,sampleMesh);

% [erroeval]=Test5_Optimization(E_ini, thickness_ini, deleteFile, sampleDim, sampleMesh);

fun = @(E_ini)Test5_Optimization(E_ini, thickness_ini, deleteFile, sampleDim, sampleMesh);

bestdat = fminsearch(fun,E_ini1)

%%%
%%
% This sectiton is the Second test of optimizaition
clear; close all; clc;

maxNumberIterations=100; %Maximum number of optimization iterations
maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=1e-6; %Tolerance on objective function value
parameterTolerance=1e-6; %Tolerance on parameter variation
displayTypeIterations='iter';
OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
    'MaxIter',maxNumberIterations,...
    'TolFun',functionTolerance,...
    'TolX',parameterTolerance,...
    'Display',displayTypeIterations,...
    'FinDiffRelStep',1e-2,...
    'DiffMaxChange',0.5);


deleteFile = 1;
sampleDim = [1.5 1.5 1.5]; %x y z
sampleMesh = [10 10 10 1]; %x y z_big_layers z_per_big_layer
E_ini1 = [1.174 3.522 2.348]/10; % E
thickness_ini = [0 0.675 1.5]; % thickness

ini_guess(:,1) = [1.174 3.522 2.348]; % E
ini_guess(:,2) = [0 0.675 1.5]; % thickness

global t;
global t2;
global er;
global mat;
mat = [];
er = [];
t = 0;
t2 = [];


% [X,Y]=Test5_Optimization(ini_guess,deleteFile,sampleDim,sampleMesh);

% [erroeval]=Test5_Optimization(E_ini, thickness_ini, deleteFile, sampleDim, sampleMesh);

fun = @(E_ini)Test6_Optimization(E_ini, thickness_ini, deleteFile, sampleDim, sampleMesh);

bestdat = fminsearch(fun,E_ini1,OPT_options)



%%
% This sectiton is the Third test of optimizaition
clear; close all; clc;

maxNumberIterations=100; %Maximum number of optimization iterations
maxNumberFunctionEvaluations=maxNumberIterations*10; %Maximum number of function evaluations, N.B. multiple evaluations are used per iteration
functionTolerance=1e-6; %Tolerance on objective function value
parameterTolerance=1e-6; %Tolerance on parameter variation
displayTypeIterations='iter';
OPT_options=optimset('fminsearch'); % 'Nelder-Mead simplex direct search'
OPT_options = optimset(OPT_options,'MaxFunEvals',maxNumberFunctionEvaluations,...
    'MaxIter',maxNumberIterations,...
    'TolFun',functionTolerance,...
    'TolX',parameterTolerance,...
    'Display',displayTypeIterations,...
    'FinDiffRelStep',1e-2,...
    'DiffMaxChange',0.5);


deleteFile = 1;
sampleDim = [1.5 1.5 1.5]; %x y z
sampleMesh = [10 10 10 1]; %x y z_big_layers z_per_big_layer
E_ini1 = [1.174 3.522 2.348]; % E
thickness_ini = [0 0.675 1.5]; % thickness

ini_guess(:,1) = [1.174 3.522 2.348]; % E
ini_guess(:,2) = [0 0.675 1.5]; % thickness

global t;
global t2;
global er;
global mat;
mat = [];
er = [];
t = 0;
t2 = [];

BCtype= 1; %1- symmetry, 2- face fixed Use onl1 1 if ratio =/= 0
strain = 0.2;
ratio = 0.75; %between X and Y dy/dx
rho_mat = 1.12e-9;

% [X,Y]=Test5_Optimization(ini_guess,deleteFile,sampleDim,sampleMesh);

% [erroeval]=Test5_Optimization(E_ini, thickness_ini, deleteFile, sampleDim, sampleMesh);

fun = @(E_ini)Test7_Optimization(E_ini, thickness_ini, deleteFile,sampleDim,sampleMesh, BCtype, strain, ratio, rho_mat);


bestdat = fminsearch(fun,E_ini1,OPT_options)


%%
%Ignore 

% % % % % %     imported = csvread('try1.csv');
% % % % % %     stretch_imported_old = imported(:,1);
% % % % % %     stress_imported_old = imported(:,2);
% % % % % %     
% % % % % %     stretch_imported = linspace(stretch_imported_old(1), stretch_imported_old(end), 20);
% % % % % %     stress_imported = interp1(stretch_imported_old, stress_imported_old, stretch_imported, 'pchip');
% % % % % %     
% % % % % %     plot(stretch_imported, stress_imported,'--k');


%%
%This section is used to compare the 2 different Possible Boundary Conditions

clear; close all; clc;
deleteFile = 1;
sampleDim = [15 15 1.5]; %x y z
sampleMesh = [20 20 10 1]; %x y z_big_layers z_per_big_layer
% E_ini = [1.174 3.522 2.348]; % E
% thickness_ini = [0 0.675 1.5]; % thickness

E_ini = [1.174 1.174 3.522 3.522 2.348 2.348]; % E
thickness_ini = [0 0.3 0.3+0.005 0.3+0.005+0.75 0.3+.0005+0.75+0.005 1.5]; % thickness

BCtype= 1; %1- symmetry, 2- face fixed Use onl1 1 if ratio =/= 0
strain = 0.2;
ratio = 0.75; %between X and Y dy/dx
rho_mat = 1.12e-9;
[X2,Y2]=UniaxialTestFunction(deleteFile,sampleDim,sampleMesh, E_ini, thickness_ini, BCtype,strain,ratio, rho_mat);

% [X1,Y1]=UniaxialTestFunction(deleteFile,sampleDim,sampleMesh, E_ini, thickness_ini, 1, strain);
% [X2,Y2]=UniaxialTestFunction(deleteFile,sampleDim,sampleMesh, E_ini, thickness_ini, 2, strain);
% 
% clf
% subplot(2,1,1);
% plot(X1,Y1);
% grid on
% grid minor
% ylabel('Stress');
% xlabel('Stretch');
% hold on
% plot(X2,Y2);
% legend('Movable','Fixed Face')
% 
% subplot(2,1,2);
% plot(X1,(Y1-Y2)*100/Y1);
% grid on
% grid minor
% ylabel('error%');
% xlabel('Stretch');

%%
E = 2.348e6;
nu = 0.45;

C10 = E/(4*(1+nu))
D1 = 2*(1+nu)*(1-2*nu)/(E*nu)