function [dydt, nonState, control] = sys_parPTO(t,y,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_parPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 11/2/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a parallel-type wave energy PTO 
% with long high- and low-pressure pipelines and an option for a throttling
% valve for enhanced pressure ripple filtering.
%
% FILE DEPENDENCY:
% ../Parallel-type PTO/
%   stateIndex_refPTO.m
% ../WEC model/
%   flapModel.m
% ../Components/
%   areaFracPWM.m
%   capAccum.m
%   deadVCap.m
%   flowCV.m
%   flowPRV.m
%
% UPDATES:
% 11/2/2023 - created from sys_refPTO.m.
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
% Define indicies of state vector
iyp_a = [];
iyp_b = [];

iyp_lin = [];
iyp_lout = [];
iyp_hin = [];
iyp_hout = [];
iyp_ro = [];

iyp_filt = [];
iy_errInt_p_filt = [];
iycontrol = [];

iytheta = [];
iytheta_dot = [];
iyrad = [];

iyLPPL = [];
iyqLP = [];
iypLP = [];
iyHPPL = [];
iyqHP = [];
iypHP = [];

stateIndex_parPTO

iWEC = [iytheta iytheta_dot iyrad];


% Calculate control input and the change in controller states (if any)
[dydt_control, control] = controller(t,y,par);

% Calculate non-state variables like coeffs., forces, and flow rates
nonState = nonStateVars(t,y,control,par);

% Calculate the hydrodynamics WEC state derivatives and output nonstate
% variables like forces and wave elevation
[dydt_WEC, nonState.torqueWEC, nonState.waveElev] = ...
                    flapModel(t,y(iWEC),nonState.T_pto,par);

% State derivatives
dydt = zeros(9 + par.WEC.ny_rad,1);

dydt(iyp_a) = 1/nonState.C_a*(par.D_WEC*y(iytheta_dot) - nonState.q_sv ...
                + nonState.q_ain - nonState.q_aout);
dydt(iyp_b) = 1/nonState.C_b*(-par.D_WEC*y(iytheta_dot) + nonState.q_sv ...
                + nonState.q_bin - nonState.q_bout);

dydt(iyp_l) = 1/nonState.C_l*(nonState.q_c ...
                - nonState.q_ain - nonState.q_bin ...
                + nonState.q_pm - nonState.q_ERUfeed ...
                - nonState.q_lPRV + nonState.q_hPRV + nonState.q_roPRV);
dydt(iyp_h) = 1/nonState.C_h*(nonState.q_aout + nonState.q_bout ...
                - nonState.q_pm ...
                + (~par.ERUconfig.outlet)*nonState.q_ERUfeed ...
                - nonState.q_rv - nonState.q_hPRV);
dydt(iyp_ro) = 1/nonState.C_ro*(nonState.q_rv - nonState.q_feed ...
                + par.ERUconfig.outlet*nonState.q_ERUfeed ...
                - nonState.q_roPRV);

dydt(iycontrol) = dydt_control;

dydt(iytheta) = dydt_WEC(1); % angular velocity
dydt(iytheta_dot) = dydt_WEC(2); % angular acceleration
dydt(iyrad) = dydt_WEC(3:end); % radiation damping states for WEC model

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [dydt_control, control] = controller(t,y,par)
        %% PI control of P_Hout using w_pm
         % Error
        err_p = y(iyp_filt) - par.control.p_ro_nom;
         % Feedforward
        w_ff = 0*par.control.w_pm_ctrl.min;
         % Control signal
        w_pm_nom = w_ff ...
            + (par.control.w_pm_ctrl.kp*err_p ...
            + par.control.w_pm_ctrl.ki*y(iy_errInt_p_filt));
        nomAboveMax = w_pm_nom > par.control.w_pm_ctrl.max;
        nomBelowMin = w_pm_nom < par.control.w_pm_ctrl.min;
        p_hAbovep_ro = par.rvConfig.included*(y(iyp_h) > y(iyp_ro)) ...
                    + ~par.rvConfig.included;
        control.w_pm = p_hAbovep_ro*(w_pm_nom ...
                  + nomAboveMax*(par.control.w_pm_ctrl.max - w_pm_nom) ...
                  + nomBelowMin*(par.control.w_pm_ctrl.min - w_pm_nom)) ...
                  + ~p_hAbovep_ro*par.control.w_pm_ctrl.min;

        %% Feedforward control of active RO inlet valve
        % determine ideal valve coefficient to satisfy dpdt limit
        C_ro = capAccum(y(iyp_ro),par.pc_ro,par.Vc_ro,par.f,par);
        q_perm = par.Sro*par.Aperm*(y(iyp_ro) - par.p_perm - par.p_osm);
        ERUfeedFlow = par.ERUconfig.present*par.ERUconfig.outlet;
        q_ro = q_perm*((~ERUfeedFlow)/par.Y ...
                      + ERUfeedFlow*(1-par.eta_ERUv^2*(1-par.Y))/par.Y);
        dp = y(iyp_h) - y(iyp_ro);
        control.kv_ideal = (sign(dp)*C_ro*par.control.dpdt_ROmax + q_ro) ...
                /sqrt(abs(dp));
        % satisfy conditions related to limitations, control objectives, 
        % and configuration
        p_roAbovepc_ro = y(iyp_ro) > par.pc_ro;
        rvOpen = ~p_roAbovepc_ro | ~par.rvConfig.active; % cond.s for valve full open
        control.kv_rv = ~rvOpen*max(0,min(par.kv_rv,control.kv_ideal)) ...
                       + rvOpen*par.kv_rv;
                      
        %% deriviatives for filtered signal and error integrals (w/
         % anti-wind-up)
        dydt_control = [    % filtered signal
                        (y(iyp_ro) - y(iyp_filt))...
                        /par.control.tau_pfilt;

                            % error integral for pressure control
                        (y(iy_errInt_p_filt) < 1e12 ...
                        & y(iy_errInt_p_filt) > -1e12) ...
                        *err_p];
                        
    end

    function nonState = nonStateVars(t,y,control,par)
        % Accumulator capacitance
        nonState.C_l = capAccum(y(iyp_l),par.pc_l,par.Vc_l,par.f,par);
        nonState.C_h = capAccum(y(iyp_h),par.pc_h,par.Vc_h,par.f,par);
        nonState.C_ro = capAccum(y(iyp_ro),par.pc_ro,par.Vc_ro,par.f,par);

        % WEC-driven pump
         % pumping chamber capacitance
        V_a = par.V_wecDead + par.D_WEC*(par.theta_max - y(iytheta));
        nonState.C_a = deadVCap(y(iyp_a),V_a,par);
        V_b = par.V_wecDead + par.D_WEC*(par.theta_max + y(iytheta));
        nonState.C_b = deadVCap(y(iyp_b),V_b,par);

         % Switching valve flow
        dp = y(iyp_a) - y(iyp_b);
        lin_lowdp = (abs(dp) < par.dp_svlin);
        nonState.q_sv = areaFracPWM(t,par.duty_sv,par.T_sv,par.tr_sv) ...
                        *(lin_lowdp*par.kv_svlin*dp ...
                        + ~lin_lowdp*sign(dp)*par.kv_sv*sqrt(abs(dp)));
         % check valves
        nonState.q_ain = flowCV(y(iyp_l) - y(iyp_a), ...
                         par.kvWECin,par.pc_WECin,par.dp_WECin);
        nonState.q_aout = flowCV(y(iyp_a) - y(iyp_h), ...
                         par.kvWECout,par.pc_WECout,par.dp_WECout);
        nonState.q_bin = flowCV(y(iyp_l) - y(iyp_b), ...
                         par.kvWECin,par.pc_WECin,par.dp_WECin);
        nonState.q_bout = flowCV(y(iyp_b) - y(iyp_h), ...
                         par.kvWECout,par.pc_WECout,par.dp_WECout);
        
         % Reaction torque on WEC
        delta_p_wp = y(iyp_b)-y(iyp_a);
        WECpumpPumping = y(iytheta_dot)*delta_p_wp < 0;
        WECpumpMotoring = ~WECpumpPumping;

        nonState.T_pto = par.D_WEC*delta_p_wp...
                    * (WECpumpPumping/par.eta_m_WEC ...
                    + WECpumpMotoring*par.eta_m_WEC);

        % House power pump/motor
        delta_p_pm = y(iyp_l) - y(iyp_h);
        nonState.pmPumping = control.w_pm*(delta_p_pm) >= 0;
        nonState.pmMotoring = control.w_pm*(delta_p_pm) < 0;

        volLoss_pm = par.APP.C_s*(delta_p_pm/(par.mu*abs(control.w_pm))) ...
                           + (delta_p_pm/par.beta)*(par.APP.V_r + 1);
        nonState.eta_v_pm = nonState.pmPumping*(1 - volLoss_pm) ...
                          + nonState.pmMotoring/(1 + volLoss_pm);
        nonState.q_pm = (nonState.pmPumping*nonState.eta_v_pm ...
                        + nonState.pmMotoring/nonState.eta_v_pm) ...
                        *par.D_pm*control.w_pm;

        mechLoss_pm = par.APP.C_v*par.mu*abs(control.w_pm)/delta_p_pm + par.APP.C_f;
        nonState.eta_m_pm = nonState.pmPumping/(1 + mechLoss_pm) ...
                          + nonState.pmMotoring*(1 - mechLoss_pm);
        nonState.Tpm = (nonState.pmPumping/nonState.eta_m_pm ...
                        + nonState.pmMotoring*nonState.eta_m_pm)...
                        *par.D_pm*delta_p_pm;

        % Reverse osmosis module
        nonState.q_perm = par.Sro*par.Aperm*(y(iyp_ro) - par.p_perm - par.p_osm);
        nonState.q_feed = nonState.q_perm/par.Y;

         % RO inlet valve
        dp = y(iyp_h) - y(iyp_ro);
        nonState.q_rv = control.kv_rv*( ...
                        par.rvConfig.included*sqrt(abs(dp))*sign(dp) ...
                        + ~par.rvConfig.included*dp);

         % ERU
        nonState.q_brine = nonState.q_feed - nonState.q_perm;
        nonState.q_ERUfeed = par.ERUconfig.present ...
                             *(par.eta_ERUv)^2*nonState.q_brine;

        % Charge Pump
        dP_SO = (y(iyp_l) - par.p_o) - par.cn*par.w_c^2; % difference between shut-off pressure and current pressure differential
        nonState.q_c = (dP_SO < 0)* sqrt(dP_SO/par.cq);

        % Pressure relief valves
        nonState.q_lPRV = flowPRV(y(iyp_l),par.lPRV.p_crack,par.lPRV.C);
        nonState.q_hPRV = flowPRV(y(iyp_h),par.hPRV.p_crack,par.hPRV.C);
        nonState.q_roPRV = flowPRV(y(iyp_ro),par.roPRV.p_crack,par.roPRV.C);

    end

end
