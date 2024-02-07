load("Data\parPTO_woPL\data_parPTO_accum_woPL_woRV_20240127_02.mat")
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
supTitleFontSize = 10;
subTitleFontSize = 9;
axFontSize = 8;
legFontSize = 7;
lineWidth = 1;


%% dpdt
width = 4.8; % one column: 3+9/16, two column: 7.5
height = 1.75;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

titleString = ['Maximum Rate of Change in Feed Pressure: Baseline PTO'];

p(1) = semilogy(Vc_h,1e-3*dpdt_max,'-k', ...
    'LineWidth',1);
p(1).HandleVisibility='off';
hold on

p(2) = plot(Vc_h([1 end]),1e-3*70e3*[1 1],'--k');

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('rate (kPa/s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')


title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

leg = legend('target limit', ...
             'Interpreter','latex');
leg.FontSize = legFontSize;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 30])
yLim = ylim;
ylim([ yLim(1) 1e3])

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;
% pos = ax.Position;
% adj = 0.125;
% ax.Position = [pos(1), pos(2)+adj, pos(3), pos(4)-adj];

%% Power loss
width = 4.8; % one column: 3+9/16, two column: 7.5
height = 1.75;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

titleString = ['Mean Power Loss Normalized by Mean Power Capture: Baseline PTO'];

p(1) = semilogy(Vc_h,100*(PP_genLoss +  PP_pmLoss)./PP_WEC,'Color',maroon,'lineStyle','--', ...
    'LineWidth',lineWidth);
hold on

% p(2) = semilogy(Vc_h,100*PP_hinPRV./PP_WEC,'Color',gold,'lineStyle',linestyles(3), ...
%     'LineWidth',lineWidth);
p(3) = semilogy(Vc_h,100*PP_roPRV./PP_WEC,'Color',blue,'lineStyle','-.', ...
    'LineWidth',lineWidth);
p(4) = semilogy(Vc_h,100*(PP_hinPRV + PP_roPRV + PP_genLoss +  PP_pmLoss)./PP_WEC,'-k', ...
    'LineWidth',lineWidth);

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('power loss ($\%$)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

title(titleString,...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

% leg = legend('hyd. motor & generator','PRV at WEC-driven pump','PRV at RO inlet','combined','Interpreter','latex');
leg = legend('hyd. motor & generator','PRV at RO inlet','combined', ...
             'Interpreter','latex');
leg.FontSize = legFontSize;
leg.FontName = 'Times';
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 30])
yLim = ylim;
ylim([0 yLim(2)])

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;
% pos = ax.Position;
% adj = 0.125;
% ax.Position = [pos(1), pos(2)+adj, pos(3), pos(4)-adj];