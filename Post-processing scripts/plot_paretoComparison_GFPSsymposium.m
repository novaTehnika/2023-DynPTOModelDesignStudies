load(['Data',filesep,'parPTO_woPL',filesep,'allData_accum_woPL_wPassiveRV_20240123_02.mat'])
SS = 2;
Vc_h_wPassiveRV = Vtotal_mesh(iiPareto);
L_wPassiveRV = PP_metric(iiPareto);
q_permMean_wPassiveRV = studyData(SS).q_permMean(iiPareto);

load(['Data',filesep,'parPTO_woPL',filesep,'allData_accumNoDpDt_woPL_wPassive_20240122.mat'])
SS = 2;
Vc_h_wPassiveRV_noDpDt = Vtotal_mesh(iiPareto);
L_wPassiveRV_noDpDt = PP_metric(iiPareto);
q_permMean_wPassiveRV_noDpDt = studyData(SS).q_permMean(iiPareto);

load(['Data',filesep,'parPTO_woPL',filesep,'allData_accum_woPL_wActiveRV_20240122_2.mat'])
SS = 2;
Vc_h_wActiveRV = Vtotal_mesh(iiPareto);
L_wActiveRV = PP_metric(iiPareto);
q_permMean_wActiveRV = studyData(SS).q_permMean(iiPareto);

load(['Data',filesep,'seriesPTO_woPL',filesep,'allData_data_seriesPTO_accum_woPL_20240122.mat'])
SS = 2;
Vc_h_series = Vtotal_mesh(iiPareto);
L_series = PP_metric(iiPareto);
q_permMean_series = studyData(SS).q_permMean(iiPareto);

load(['Data',filesep,'parPTO_woPL',filesep,'data_parPTO_accum_woPL_woRV_20240127_02.mat'])
SS = 2;
Vc_h_baseline = Vc_h;
L_baseline = 100*(PP_hinPRV + PP_roPRV + PP_genLoss +  PP_pmLoss)./PP_WEC;
q_permMean_baseline = q_permMean;

Vc_h_woRV = 17.6;
L_limit = 5;

SS = 2;

black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
pink = [255 192 203]/256;
blue1 = [0 150 255]/256;
purple = [128 0 128]/256;
green1 = [0 255 150]/256;
color = [maroon; gold; blue; orange; green; pink; blue1; purple; green1];

linestyles = {'-', '--', '-.', ':','-', '--', '-.', ':',};

bottomEdge = 1;
leftEdge = 3;
width = 4.8; % one column: 3+9/16, two column: 7.5
height = 2.5;
supTitleFontSize = 9;
subTitleFontSize = 9;
axFontSize = 8;
legFontSize = 8;
lineWidth = 0.5;
markerSize = 4;


clearvars leg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power loss

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ip = 1;
p(ip) = semilogy(Vc_h_woRV*[1 1],[0.01 100],'-.k','LineWidth',lineWidth);
hold on
ip = ip+1;
p(ip) = semilogy(Vc_h_wPassiveRV,L_wPassiveRV,'Color',maroon,'Marker','x',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
ip = ip+1;
p(ip) = semilogy(Vc_h_wActiveRV,L_wActiveRV,'Color',gold,'Marker','o',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
ip = ip+1;
p(ip) = semilogy(Vc_h_series,L_series,'Color',blue,'Marker','^',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
ip = ip+1;
p(ip) = semilogy(Vc_h_baseline,L_baseline,'Color',black,...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
ip = ip+1;
p(ip) = semilogy(Vc_h_wPassiveRV_noDpDt,L_wPassiveRV_noDpDt,'Color',green,'Marker','^',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
ip = ip+1;
p(ip) = semilogy(xlim,L_limit*[1 1],'--k',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);


grid on

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Mean Power Loss Normalized to Mean Power Capture',newline, ...
                'vs. Total High-Pressure Accumulator Volume'];
title(titleString,...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

ylim([1 100])
xlim([0 20])

leg = legend('baseline','passive elem.','active elem.', ...
    'series-type','baseline, no constraint','passive elem., no constraint','5% loss limit', ...
             'Interpreter','latex');
leg.FontSize = legFontSize;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Permeate Production

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ip = 1;
p(ip) = plot(Vc_h_wPassiveRV,24*3600*q_permMean_wPassiveRV,'Color',maroon,'Marker','x',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
hold on
ip = ip+1;
p(ip) = plot(Vc_h_wActiveRV,24*3600*q_permMean_wActiveRV,'Color',gold,'Marker','o',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
ip = ip+1;
p(ip) = plot(Vc_h_series,24*3600*q_permMean_series,'Color',blue,'Marker','^',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
ip = ip+1;
p(ip) = plot(Vc_h_baseline,24*3600*q_permMean_baseline,'Color',black,...
    'LineWidth',lineWidth,'MarkerSize',markerSize);
ip = ip+1;
p(ip) = plot(Vc_h_wPassiveRV_noDpDt,24*3600*q_permMean_wPassiveRV_noDpDt,'Color',green,'Marker','^',...
    'LineWidth',lineWidth,'MarkerSize',markerSize);


grid on

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('permeate rate (m$^3$/day)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Mean Permeate Production',newline, ...
                'vs. Total High-Pressure Accumulator Volume'];
title(titleString,...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

% ylim([1 100])
xlim([0 20])

leg = legend('passive elem.','active elem.', ...
    'series-type','baseline, no constraint','passive elem., no constraint', ...
             'Interpreter','latex');
leg.FontSize = legFontSize;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')