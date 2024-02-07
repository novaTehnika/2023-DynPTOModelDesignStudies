X = Vc_h;
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
height = 7;
fontSize = 9;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 3;

iax = 1;
ax(iax) = subplot(n_plots,1,1);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;



p(iax,1) = semilogy(X,Y,'-k', ...
    'LineWidth',1);
p(iax,1).HandleVisibility='off';
hold on

p(iax,2) = plot(X([1 end]),1e-3*par.control.dpdt_ROmax*[1 1],'--k')
legLabels(1) = convertCharsToStrings( ...
        ['target limit']);

xlabel('volume (1000L)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('rate of change in pressure (kPa/s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title([varTitle,' Without Ripple Control Valve: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg(iax) = legend(legLabels);
leg(iax).FontSize = fontSize-1;
leg(iax).FontName = 'Times';
% rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg(iax), 'Location', 'best')
% set(leg, 'Location', 'southoutside')
% xlim([0 max(Vtotal)])
yLim = ylim;
ylim([0 yLim(2)])
% ylim([0 1e-3*1.5*par.control.dpdt_ROmax])

% %%%%%%%%%%%% power loss
Y = 100*(PP_hinPRV + PP_roPRV + PP_genLoss +  PP_pmLoss)./PP_WEC;
iax = 2;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;


legLabels(1) = convertCharsToStrings( ...
        ['combined']);
p(iax,1) = semilogy(X,Y,'-k', ...
    'LineWidth',1);
hold on

legLabels(2) = convertCharsToStrings( ...
        ['WEC outlet PRV']);
p(iax,2) = semilogy(X,100*PP_hinPRV./PP_WEC,'-r', ...
    'LineWidth',1);

legLabels(3) = convertCharsToStrings( ...
        ['RO feed PRV']);
p(iax,3) = semilogy(X,100*PP_roPRV./PP_WEC,'-g', ...
    'LineWidth',1);

legLabels(4) = convertCharsToStrings( ...
        ['pump/motor and generator']);
p(iax,4) = semilogy(X,100*(PP_genLoss +  PP_pmLoss)./PP_WEC,'k','Marker','x', ...
    'LineWidth',1);


xlabel('volume (1000L)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power loss ($\%$)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Portion of WEC Power lost to PRVs: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg(iax) = legend(legLabels);
leg(iax).FontSize = fontSize-1;
leg(iax).FontName = 'Times';
% rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg(iax), 'Location', 'best')
% set(leg, 'Location', 'southoutside')
% xlim([0 max(Vtotal)])
yLim = ylim;
ylim([0 yLim(2)])
% ylim([0 1e-3*1.5*par.control.dpdt_ROmax])


% %%%%%%%%%%%% permeate production
Y = q_permMean*3600*24;
iax = 3;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;


p(iax,1) = plot(X,Y,'-k', ...
    'LineWidth',1);

xlabel('volume (1000L)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('permeate production ($m^3/day$)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Portion of WEC Power lost to PRVs: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

% leg = legend(legLabels);
% leg.FontSize = fontSize-1;
% leg.FontName = 'Times';
% % rect = [0.5, -0.2, 0.25, 0.15];
% % set(leg, 'Position', rect)
% set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
% xlim([0 max(Vtotal)])
yLim = ylim;
% ylim([0 yLim(2)])
% ylim([0 1e-3*1.5*par.control.dpdt_ROmax])

linkaxes(ax,'x')