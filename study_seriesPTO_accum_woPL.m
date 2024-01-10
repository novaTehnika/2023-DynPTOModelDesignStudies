% study_seriesPTO_accum_woPL.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 1/9/2024
%
% PURPOSE/DESCRIPTION:
% This script performs parameter variation studies
% using the model contained in sys_seriesPTO.m and solved by
% sim_seriesPTO.m.
% The parameter initiallization functions are called within this
% script before the sim_seriesPTO.m script is called.
%
% This specific script studies the size of high-pressure accumulators and
% as they effect pressure rate of change at the RO module and the power 
% losses in the system (e.g. power loss through pressure relief valves.
%
% This script is set up to be run as part of a SLURM job array. The
% following lines are required before this script is called:
%   iVar = ${SLURM_ARRAY_TASK_ID};
%   SS=1;
%
% FILE DEPENDENCY:
% ./Series-type PTO/
%   initialConditionDefault_seriesPTO
%   parameters_seriesPTO.m
%   sim_seriesPTO.m
%   stateIndex_seriesPTO.m
%   sys_seriesPTO.m
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
% ./Components/Pipeline
%   flowR.m
%   lineCap.m
%   pipelineNPi.m
% ./Utilities/
%   startParPool.m
%   statsTimeVar_cdf.m
%   get_current_git_hash.m
%   leadingZeros.m
%
% UPDATES:
% 1/9/2024 - Created from study_parPTO_accum_wPassiveRV.m
%
% Copyright (C) 2024  Jeremy W. Simmons II
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
% clear
% clc
addpath('WEC model') 
addpath(['WEC model' filesep 'WECdata']) 
addpath('Series-type PTO')
addpath('Components')
addpath(['Components' filesep 'Pipeline'])
addpath('Sea States')
addpath('Solvers')
addpath('Utilities')
[git_hash_string, git_status_string] = get_current_git_hash();

%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.tstart = 0; %[s] start time of simulation
par.tend = 1000;%2000; %[s] end time of simulation

par.Tramp = 250; % [s] excitation force ramp period
par.TrampWEC = min(25,par.Tramp); % [s] excitation force ramp period

% Solver parameters
par.solver = 'fixed time'; % 'variable time' OR 'fixed time'
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
par = parameters_seriesPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

%% Special modifications to base parameters
% par.Sro = 3700; % [m^3]
% par.D_WEC = 0.23;         % [m^3/rad] flap pump displacement

% Operating parameters
p_h_nom = 1e6*[999, 7.4049, 10.3516, 7.6581, 10.3516,9.3226]; % [Pa]
p_ro_nom = 1e6*[999, 5.0533, 8, 5.3066, 8, 6.9710]; % [Pa]
w_c = [2500 2500 3000 3000 3000 3000]*2*pi/60; % [(rpm) -> rad/s]

par.control.p_h_nom = p_h_nom(SS);
par.control.p_ro_nom = p_ro_nom(SS);
par.w_c = w_c(SS);
par.duty_sv = 0;

% Configuration
par.ERUconfig.present = 1;

par.plConfig.included = 0;

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total accumulator volume
% distribution between motor inlet and RO inlet

nVar1 = 5;
Vtotal = 1e-3*logspace(log10(2e3),log10(30e3),nVar1);% [L->m^3] total accumulator volume
nVar2 = 9;
X = linspace(0.1,0.5,nVar2); % [-] accumulator volume distribution 1 - all at RO inlet, 0 - all at motor inlet
nVar3 = 4;
D_pm = linspace(1000,4000,nVar3)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad]  Motor displacement

[meshVar.Vtotal, meshVar.X, meshVar.D_pm] = meshgrid(Vtotal,X,D_pm);
Vtotal_mesh = meshVar.Vtotal(:);
X_mesh = meshVar.X(:);
D_pm_mesh = meshVar.D_pm(:);

nVar = length(Vtotal_mesh);

saveSimData = 0; % save simulation data (1) or just output variables (0)

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change design parameter
Vc_h = (1 - X_mesh(iVar))*Vtotal_mesh(iVar);
par.Vc_hin = 0.5*Vc_h;
par.Vc_hout = 0.5*Vc_h;
par.Vc_ro = X_mesh(iVar)*Vtotal_mesh(iVar);
par.D_pm = D_pm_mesh(iVar);

% Define state indices
par.iy = stateIndex_seriesPTO(par);

% Define initial conditions
y0 = initialConditionDefault_seriesPTO(par); % default ICs, provides 'y0'

% run simulation
ticSIM = tic;
out = sim_seriesPTO(y0,par);
toc(ticSIM)

% Calculate metrics
% max rate of change in pressure
% 97th percentile rate of change
% power loss from valve
% power loss through PRVs
% permeate production
% power loss from pump/motor and power generated for normalization
% 
q_permMean = mean(out.q_perm);
PP_WEC = mean(out.power.P_WEC);
PP_wp = mean(out.power.P_wp);
PP_pmLoss = mean(out.power.P_pmLoss);
PP_gen = mean(out.power.P_gen);
PP_hinPRV = mean(out.power.P_hinPRV);
PP_roPRV = mean(out.power.P_roPRV);
dpdt_max = max(abs(out.dydt(:,par.iy.p_ro)));

try
    dist_dpdt = statsTimeVar_cdf(out.t,abs(out.dydt(:,par.iy.p_ro)));
    dpdt_97(iVar) = dist_dpdt.xi(find(dist_dpdt.f > 0.97,1,'first'));
catch
    dist_dpdt = nan;
    dpdt_97(iVar) = nan;
end

if ~saveSimData
    clear out
end

%% %%%%%%%%%%%%   End Computations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poolobj = gcp('nocreate'); delete(poolobj);

%% %%%%%%%%%%%%   Save Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeStamp = datetime("now",'format','yyyy-MM-dd''T''HH:mm'); % time in ISO8601

% Save data
filename = ['data_seriesPTO_accum_woPL', ...
            '_',char(datetime("now",'Format','yyyyMMdd')), ...
            '_',num2str(SS,leadingZeros(999)), ...
            '_',num2str(iVar,leadingZeros(nVar))];
save(filename,'-v7.3')

return
