X = 1e3*Vc_h;
par = parBase;
% select variable to plot
maxOr97 = 1;
switch maxOr97
  case 1
    Y = 1e-3*dpdt_max;
    varTitle = 'Maximum Rate of Change in Pressure';
  case 2
    Y = 1e-3*dpdt_97;
    varTitle = '97th Percentile Rate of Change in Pressure';
end

black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];

linestyles = {'-', '--', ':', '-.'};

bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 3.25;
fontSize = 9;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 1;
ax1 = subplot(n_plots,1,1);
ax1.FontName = 'times';
ax1.FontSize = fontSize-1;

hold on

plot(X([1 end]),1e-3*par.control.dpdt_ROmax*[1 1],'--k')
legLabels(1) = convertCharsToStrings( ...
        ['target limit']);

p = plot(X,Y,'-k', ...
    'LineWidth',1);

xlabel('volume (L)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('rate of change in pressure (kPa/s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title([varTitle,' Without Ripple Control Valve: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels);
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
% rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
% xlim([0 max(Vtotal)])
yLim = ylim;
ylim([0 yLim(2)])
% ylim([0 1e-3*1.5*par.control.dpdt_ROmax])
