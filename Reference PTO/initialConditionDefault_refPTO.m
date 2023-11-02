% initialConditionDefault_refPTO.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/12/2023
%
% PURPOSE/DESCRIPTION:
% This script loads default intial conditions for the refPTO model.
%
% FILE DEPENDENCY:
%
% UPDATES:
% 6/12/2023 - Created from initialConditionDefault_parallelPTO.m.
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

y0(1,iytheta) = 0.1*pi/2; % [rad]
y0(1,iytheta_dot) = 0; % [rad/s]
y0(1,iyrad) = zeros(1,length(iyrad));

y0(1,iyp_a) = 1e6; % [Pa]
y0(1,iyp_b) = 1e6; % [Pa]

y0(1,iyp_l) = par.control.p_l_nom; % [Pa]
y0(1,iyp_h) = par.control.p_h_nom; % [Pa]
y0(1,iyp_ro) = par.control.p_h_nom; % [Pa]

y0(1,iyp_filt) = par.control.p_h_nom;
y0(1,iy_errInt_p_filt) = 0;

