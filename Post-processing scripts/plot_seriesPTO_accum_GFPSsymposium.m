load("Data\seriesPTO_woPL\allData_data_seriesPTO_accum_woPL_20240122.mat")
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
width = 4.5; % one column: 3+9/16, two column: 7.5
height = 3.75;
supTitleFontSize = 10;
subTitleFontSize = 9;
axFontSize = 8;
legFontSize = 7;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

titleString = ['Series-Type Archetecture Performance',newline,...
                'as a Function of Total Accumulator Volume'];
sgtitle(titleString,...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

n_plots = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power loss
iax = 1;
ax(iax) = subplot(n_plots,1,1:2);

ip = 1;
p(ip,iax) = semilogy(V_metric_opt,100*(studyData(SS).L_pmLoss(iiPareto)+studyData(SS).L_genLoss(iiPareto)),'Color',maroon,'Marker','o');
hold on
ip = ip+1;
p(ip,iax) = semilogy(V_metric_opt,100*(studyData(SS).L_hinPRV(iiPareto)),'Color',gold,'Marker','^');
ip = ip+1;
p(ip,iax) = semilogy(V_metric_opt,100*(studyData(SS).L_roPRV(iiPareto)),'Color',blue,'Marker','square');
ip = ip+1;
p(ip,iax) = semilogy(V_metric_opt,PP_metric(iiPareto),'k','Marker','x');

grid on

ax(iax).FontName = 'times';
ax(iax).FontSize = axFontSize;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Mean Power Loss Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

leg = legend('hyd. motor & generator','PRV at WEC-driven pump','PRV at RO inlet','combined');
leg.FontSize = legFontSize;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% distibution of volume
iax = 2;
ax(iax) = subplot(n_plots,1,3);

hold on

ip = 1;
p(ip,iax) = plot(V_metric_opt,X_opt,'k','Marker','x');
ip = ip+1;

grid on

ax(iax).FontName = 'times';
ax(iax).FontSize = axFontSize;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('portion', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Portion of Volume at RO Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

linkaxes(ax,'x')