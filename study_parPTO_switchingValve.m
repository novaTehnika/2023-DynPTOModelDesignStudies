% study_parPTO_switchingValve.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 11/22/2023
%
% PURPOSE/DESCRIPTION:
% This script performs parameter variation studies
% using the model contained in sys_parPTO.m and solved by
% sim_parPTO.m.
% The parameter initiallization functions are called within this
% script before the sim_parPTO.m script is called.
%
% This specific script studies the size of swithcing valve at the 
% WEC-driven pump.
%
% This script is set up to be run as part of a SLURM job array. The
% following lines are required before this script is called:
%   iVar = ${SLURM_ARRAY_TASK_ID};
%   SS=1;
%
% FILE DEPENDENCY:
% ./Parallel-type PTO/
%   initialConditionDefault_parPTO
%   parameters_parPTO.m
%   sim_parPTO.m
%   stateIndex_parPTO.m
%   sys_parPTO.m
% ./WEC model/
%   flapModel.m
%   hydroStaticTorque.m
%   parameters_WECmodel.m
% ./WEC model/WECdata
%   nemohResults_vantHoff2009_20180802.mat
%   vantHoffTFCoeff.mat
% ./Solvers/
%   deltaE_NI.m
%   deltaV_NI.m
%   ode1.m
% ./Components/
%   areaFracPWM.m
%   capAccum.m
%   deadVCap.m
%   flowCV.m
%   flowPRV.m
% ./Utilities/
%   get_current_git_hash.m
%   leadingZeros.m
%
% UPDATES:
% 11/22/2023 - Created from study_refPTO_switchingValve.m.
%
% Copyright (C) 2023  Jeremy W. Simmons II
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% clc
addpath('WEC model') 
addpath(['WEC model' filesep 'WECdata']) 
addpath('Reference PTO')
addpath('Components')
addpath('Sea States')
addpath('Solvers')
addpath('Utilities')
git_hash_string = get_current_git_hash();

%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 500; %[s] end time of simulation

% Solver parameters
% par.odeSolverRelTol = 1e-4; % Rel. error tolerance parameter for ODE solver
% par.odeSolverAbsTol = 1e-4; % Abs. error tolerance parameter for ODE solver
par.MaxStep = 5e-5;             % [s] for fixed time solver, this is the step size for solver
par.downSampledStepSize = 1e-2; % [s] specifies time step for data output
if mod(par.downSampledStepSize,par.MaxStep)
    warning('down-sampled time step is not an integer multiple of the maximum step size')
end

% Sea State and Wave construction parameters
Hs = [2.34 2.64 5.36 2.05 5.84 3.25];
Tp = [7.31 9.86 11.52 12.71 15.23 16.5];
par.wave.Hs = Hs(SS);
par.wave.Tp = Tp(SS);
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition
par.wave.rngSeedPhase = 3; % seed for the random number generator

% load parameters
par = parameters_parPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

% Define initial conditions
stateIndex_parPTO % load state indices, provides 'iy_...'
initialConditionDefault_parPTO % default ICs, provides 'y0'

%% Special modifications to base parameters
% par.Sro = 3700; % [m^3]
% par.D_WEC = 0.23;         % [m^3/rad] flap pump displacement
p_ro_nom = 1e6*[4.0000 4.9435 8.0000 5.2661 8.0000 7.1052]; % [Pa]
par.control.p_ro_nom = p_ro_nom(SS); % [Pa]

par.ERUconfig.present = 1;
par.ERUconfig.outlet = 1;

par.rvIncluded = 0; % RO inlet valve is 1 - present, 0 - absent
par.rvConfig = (0)*par.rvIncluded; % RO inlet valve is 1 - active, 0 - passive
% dp_rated = 1e5; % [Pa] 
% q_rated = (100)*60/1e3; % [(lpm) -> m^3/s]
% par.kv_rv = q_rated/dp_rated;

par.duty_sv = 0.25;

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nVar1 = 10;
kv = logspace(log10(1e-4),log10(1e-3),nVar1);% [m^3/s/Pa] valve coefficient for swithcing valve, full open
nVar2 = 10;
Tsw = logspace(log10(0.1),log10(5),nVar2);% [s] switching period

[meshVar.kv, meshVar.Tsw] = meshgrid(kv,Tsw);
kv_mesh = meshVar.kv(:);
Tsw_mesh = meshVar.Tsw(:);

nVar = length(kv_mesh);

saveSimData = 1; % save simulation data (1) or just output variables (0)

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change design parameter
par.kv_sv = kv_mesh(iVar);
par.T_sv = Tsw_mesh(iVar);

% run simulation
ticSIM = tic;
out = sim_parPTO(y0,par);
toc(ticSIM)

% Calculate metrics
it_vec = find(out.t>=par.tstart);
PP_WEC(iVar) = mean(out.power.P_WEC(it_vec));
PP_wp(iVar) = mean(out.power.P_wp(it_vec));
eff_wecPump(iVar) = PP_wp(iVar)/PP_WEC(iVar);

if saveSimData
    simOut(iVar) = out;
end

%% %%%%%%%%%%%%   Post-process data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = length(kv);
K = length(Tsw);
PP_WEC_2D = zeros(J,K);
eff_wecPump_2D = zeros(J,K);

test = 1;
for j = 1:J
    for k = 1:K
        i = K*(j-1) + k;
        PP_WEC_2D(j,k) = PP_WEC(i);
        PP_wp_2D(j,k) = PP_wp(i);
        eff_wecPump_2D(j,k) = eff_wecPump(i); 
        test = (meshVar.Tsw(i) == Tsw(k)) && (meshVar.kv(i) == kv(j)) && test;
    end
end

if ~test; error('indexing incorrect'); end
clearvars test

%% %%%%%%%%%%%%   Save Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save: time in ISO8601
filename = ['data_parPTO_switchingValve', ...
            '_',char(datetime("now",'Format','yyyyMMdd')), ...
            '_',num2str(SS,leadingZeros(999)), ...
            '_',num2str(iVar,leadingZeros(nVar))];
save(filename,'-v7.3')

%% %%%%%%%%%%%%   End Computations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WEC pump efficiency
X = kv*sqrt(1000)*1e3;
Y = Tsw;
Z = eff_wecPump_2D; 
inc = 250;
levels = round(max(eff_wecPump_2D(:)) - logspace(log10(0.01),log(0.9),8),2);

figure
C1 = contour(log10(X),Y,Z,levels,'-k','ShowText','on');
xt = xticks;
xticklabels("10^{"+string(xt')+"}");

xlabel('flow coefficient (L/s/kPa^{1/2})');
ylabel('switching period (s)');
title('WEC-driven Pump Efficiency')

% Absorbed power
X = kv*sqrt(1000)*1e3;
Y = Tsw;
Z = PP_WEC_2D*1e-3;


figure
C2 = contour(log10(X),Y,Z,'-k','ShowText','on');
xt = xticks;
xticklabels("10^{"+string(xt')+"}");

xlabel('flow coefficient (L/s/kPa^{1/2})');
ylabel('switching period (s)');
title('WEC Power Absorption (kW)')

% Pump output power
X = kv*sqrt(1000)*1e3;
Y = Tsw;
Z = PP_wp_2D*1e-3;


figure
C2 = contour(log10(X),Y,Z,'-k','ShowText','on');
xt = xticks;
xticklabels("10^{"+string(xt')+"}");

xlabel('flow coefficient (L/s/kPa^{1/2})');
ylabel('switching period (s)');
title('WEC-Driven Pump Power Output (kW)')
