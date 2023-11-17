% stateIndex_parPTO.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 11/2/2023
%
% PURPOSE/DESCRIPTION:
% This script loads the state indices for the parPTO model
%
% FILE DEPENDENCY:
%
% UPDATES:
% 11/2/2023 - Created from stateIndex_refPTO.m.
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

iyp_a = 1;
iyp_b = 2;

iyp_lin = 3;
iyp_lout = 4;
iyp_hin = 5;
iyp_hout = 6;
iyp_ro = 7;

iyp_filt = 8;
iy_errInt_p_filt = 9;
iycontrol = [iyp_filt; iy_errInt_p_filt];

iytheta = 10;
iytheta_dot = 11;
iyrad = (1:par.WEC.ny_rad) + iytheta_dot; % state vector indices for radiation damping states for WEC model

iyLPPL = (1:2*par.n_seg(1)-1) + iyrad(end);
iyqLP = (1:2:2*par.n_seg(1)-1) + (iyLPPL(1)-1);
iypLP = (2:2:2*par.n_seg(1)-1) + (iyLPPL(1)-1);

iyHPPL = (1:2*par.n_seg(2)-1) + iyLPPL(end);
iyqHP = (1:2:2*par.n_seg(2)-1) + (iyHPPL(1)-1);
iypHP = (2:2:2*par.n_seg(2)-1) + (iyHPPL(1)-1);

ny = iyHPPL(end);