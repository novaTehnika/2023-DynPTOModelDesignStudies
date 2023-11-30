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
% w_c = [3000 3000 3000 3000 3000 3000]*2*pi/60; % [(rpm) -> rad/s]
par.control.p_ro_nom = p_ro_nom(SS);
par.duty_sv = 0;

% Configuration
par.ERUconfig.present = 1;
par.ERUconfig.outlet = 1;

par.rvConfig.included = 0; % RO inlet valve is 1 - present, 0 - absent
par.rvConfig.active = (0)*par.rvConfig.included; % RO inlet valve is 1 - active, 0 - passive

par.pc_lin = 0.15e6; % [Pa] charge pressure
par.pc_lout = 0.15e6; % [Pa] charge pressure


%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total low-pressure accumulator volume
nVar1 = 20;
Vc_l = 1e-3*logspace(log10(500),log10(5000),nVar1); % [(L) -> m^3]

% portion of low-pressure accumulator volume at WEC-driven pump inlet
X = 0.1:0.1:0.9;
nVar2 = numel(X);

% Diameter of low-pressure pipeline
nVar3 = 10;
d_LPPL = logspace(log10(0.05),log10(0.4),nVar3); % [m]

% charge pump speed
w_c = (1700:100:3500)*2*pi/60; % [(rpm) -> rad/s]
nVar4 = numel(w_c);

[meshVar3D.Vc_l, meshVar3D.X, meshVar3D.d_LPPL] = meshgrid(Vc_l,X,d_LPPL);
Vc_l_mesh = meshVar3D.Vc_l(:);
X_mesh = meshVar3D.X(:);
d_LPPL_mesh = meshVar3D.d_LPPL(:);

nVar = numel(Vc_l_mesh);

saveSimData = 0; % save simulation data (1) or just output variables (0)

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% change design parameters
 % accumulator volume
par.Vc_lin = (1-X_mesh(iVar))*Vc_l_mesh(iVar);  % inlet of LP pipeline
par.Vc_lout = X_mesh(iVar)*Vc_l_mesh(iVar);     % outlet of LP pipeline
 % pipeline diameter
par.d_line(1) = d_LPPL_mesh(iVar);

parBase = par; clearvars par
parfor iw_c = 1:nVar4

    par = parBase;
    par.w_c = w_c(iw_c);

    % run simulation
    ticSIM = tic;
    out = sim_parPTO(y0,par);
    toc(ticSIM)

    % postprocess simulation data
    if sum(imag(abs(out.p_lout) ...
            + abs(out.p_a) ...
            + abs(out.p_b))) == 0

         % Mean pressure at WEC-driven pump inlet
        p_loutMean(iw_c) = mean(out.p_lout);
         % Variation in pressure at WEC-driven pump inlet
        p_loutMax(iw_c) = max(out.p_lout);
        p_loutMin(iw_c) = min(out.p_lout);
        p_loutVar(iw_c) = var(out.p_lout);
        p_loutStd(iw_c) = std(out.p_lout);
         % Minimum pressure in WEC-driven pump chambers
        p_wpMin(iw_c) = min(min(out.p_a),min(out.p_b));

         % Electric power consumption of charge pump
        P_cElec(iw_c) = mean(out.power.P_cElec);
        P_cElec_norm(iw_c) = P_cElec(iw_c)/mean(out.power.P_WEC);
         % Power losses from charge pump
        P_cLoss(iw_c) = mean(out.power.P_cLoss);
        L_c(iw_c) = P_cLoss(iw_c)/mean(out.power.P_WEC);

         % power loss from pipeline
        P_LPPL(iw_c) = out.power.P_LPPL;
        L_LPPL(iw_c) = P_LPPL(iw_c)/mean(out.power.P_WEC);

    else

        p_loutMean(iw_c) = nan;
         % Variation in pressure at WEC-driven pump inlet
        p_loutMax(iw_c) = nan;
        p_loutMin(iw_c) = nan;
        p_loutVar(iw_c) = nan;
        p_loutStd(iw_c) = nan;
         % Minimum pressure in WEC-driven pump chambers
        p_wpMin(iw_c) = nan;

         % Electric power consumption of charge pump
        P_cElec(iw_c) = nan;
        P_cElec_norm(iw_c) = nan;
         % Power losses from charge pump
        P_cLoss(iw_c) = nan;
        L_c(iw_c) = nan;

         % power loss from pipeline
        P_LPPL(iw_c) = nan;
        L_LPPL(iw_c) = nan;

    end

    if saveSimData
        outSave(iw_c) = out; %#ok<UNRCH>
    end

end
par = parBase; clearvars parBase


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
