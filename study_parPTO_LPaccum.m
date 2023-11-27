% study_parPTO_LPaccum.m script m-file
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
% This specific script studies the size of the low-pressure accumulators 
% and the size of the low-pressure pipeline.
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
% 11/27/2023 - Created from study_parPTO_chargePumpAccum.m
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
% clear
% clc
addpath('WEC model') 
addpath(['WEC model' filesep 'WECdata']) 
addpath('Parallel-type PTO')
addpath('Components')
addpath(['Components' filesep 'Pipeline'])
addpath('Sea States')
addpath('Solvers')
addpath('Utilities')
git_hash_string = get_current_git_hash();

%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.Tramp = 250; % [s] excitation force ramp period
par.tstart = 0; %[s] start time of simulation
par.tend = 2000; %[s] end time of simulation

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
w_c = [3000 3000 3000 3000 3000 3000]*2*pi/60; % [(rpm) -> rad/s]
par.control.p_ro_nom = p_ro_nom(SS);
par.duty_sv = 0;

% Configuration
par.ERUconfig.present = 1;
par.ERUconfig.outlet = 1;

par.rvConfig.included = 1; % RO inlet valve is 1 - present, 0 - absent
par.rvConfig.active = (0)*par.rvConfig.included; % RO inlet valve is 1 - active, 0 - passive
par.kv_rv = (4.7)/sqrt(1e3)/1000; % [(L/s/kPa^0.5) -> m^3/s/Pa^0.5]

par.Vc_hin = (2500)*1e-3; % [(L) -> m^3] gas volume at charge pressure
par.Vc_hout = (2500)*1e-3; % [(L) -> m^3] gas volume at charge pressure
par.Vc_ro = (5000)*1e-3; % [(L) -> m^3] gas volume at charge pressure

par.pc_l = 0.15e6; % [Pa] charge pressure

par.D_pm = (1000)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad]  Motor displacement
par.w_pm_max = (1750)/60*2*pi; % [(rpm) -> rad/s] maximum speed of motor

% charge pump pump curve parameters
par.cn = 5.5;
par.cq = -6e6;

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total low-pressure accumulator volume
nVar1 = 20;
Vc_l = 1e-3*logspace(log10(100),log10(5000),nVar1); % [(L) -> m^3]

% portion of low-pressure accumulator volume at WEC-driven pump inlet
X = 0.1:0.1:0.9;
nVar2 = numel(X);

% Diameter of low-pressure pipeline
nVar3 = 10;
d_LPPL = logspace(log10(0.05),log10(0.4),nVar3); % [m]


[meshVar.Vc_l, meshVar.X, meshVar.d_LPPL] = meshgrid(Vc_l,X,d_LPPL);
Vc_l_mesh = meshVar.Vc_l(:);
X_mesh = meshVar.X(:);
d_LPPL_mesh = meshVar.d_LPPL(:);

nVar = numel(Vc_l_mesh);

saveSimData = 1; % save simulation data (1) or just output variables (0)

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change design parameters
 % accumulator volume
par.Vc_lin = (1-X_mesh(iVar))*Vc_l_mesh(iVar);  % inlet of LP pipeline
par.Vc_lout = X_mesh(iVar)*Vc_l_mesh(iVar);     % outlet of LP pipeline
 % pipeline diameter
par.d_line(1) = d_LPPL_mesh(iVar);

% run simulation
ticSIM = tic;
out = sim_parPTO(y0,par);
toc(ticSIM)

if ~saveSimData
    clear out
end

%% %%%%%%%%%%%%   Save Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeStamp = datetime("now",'format','yyyy-MM-dd''T''HH:mm'); % time in ISO8601

% Save data
filename = ['data_parPTO_LPaccum', ...
            '_',char(datetime("now",'Format','yyyyMMdd')), ...
            '_',num2str(SS,leadingZeros(999)), ...
            '_',num2str(iVar,leadingZeros(nVar))];
save(filename,'-v7.3')

%% %%%%%%%%%%%%   End Computations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poolobj = gcp('nocreate'); delete(poolobj);

return
