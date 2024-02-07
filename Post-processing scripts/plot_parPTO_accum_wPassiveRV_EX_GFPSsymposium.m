load("Data\data_parPTO_accum_woPL_wPassiveRV_20240131_02_9090L.mat")

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

supTitleFontSize = 9;
subTitleFontSize = 9;
axFontSize = 8;
bottomEdge = 1;
leftEdge = 3;
lineWidth = 0.5;

%% Pressure and control
width = 4.8; % one column: 3+9/16, two column: 7.5
height = 1.75;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

titleString = ['Pressure and Control: Parallel-Type with Passive Element'];

set(fig,'defaultAxesColorOrder',[black; black]);

yyaxis left
plot(out.t([1 end]),1e-6*out.par.control.p_ro_nom*[1 1],'--k','LineWidth',lineWidth)
hold on
plot(out.t,1e-6*out.p_ro,'-','LineWidth',lineWidth,'Color',maroon)
plot(out.t,1e-6*out.p_hin,'-.','LineWidth',lineWidth,'Color',gold)
plot(out.t,1e-6*out.p_ro,'-','LineWidth',lineWidth,'Color',maroon,'HandleVisibility','off')


yLim = ylim;
ylim([4 yLim(2)])

xlabel('time (s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('pressure (MPa)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

yyaxis right
plot(out.t,1e3*out.q_pm(:),'LineStyle',':','LineWidth',lineWidth*2,'Color',blue)
ylabel('flow rate (L/s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

yLim = ylim;
ylim([0 yLim(2)])

leg = legend('nom. feed pressure','$p_{ro}$','$p_{h}$','$q_{m}$', ...
             'Interpreter','latex');
leg.FontSize = axFontSize;
leg.FontName = 'Times';
leg.Orientation = 'horizontal';
pos = leg.Position;
leg.Position = [0.5-pos(3)/2, 0, pos(3), pos(4)];
leg.Box = 'off';

title(titleString,...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;
pos = ax.Position;
adj = 0.125;
ax.Position = [pos(1), pos(2)+adj, pos(3), pos(4)-adj];


%% Rate of change in pressure
width = 4.8; % one column: 3+9/16, two column: 7.5
height = 1.75;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

titleString = ['Rate of Change in Feed Pressure: Parallel-Type with Passive Element'];

plot(out.t,1e-3*out.dydt(:,par.iy.p_ro),'LineWidth',lineWidth,'Color',maroon)
hold on



plot(out.t([1 end]),70*[1 1],'--k','LineWidth',lineWidth)
p(3) = plot(out.t([1 end]),70*[-1 -1],'--k','LineWidth',lineWidth);
p(3).HandleVisibility='off';

addpath('Utilities')
dist_dpdt = statsTimeVar_cdf(out.t,abs(out.dydt(:,par.iy.p_ro)));
dpdt_97 = dist_dpdt.xi(find(dist_dpdt.f > 0.97,1,'first'));
plot(out.t([1 end]),1e-3*dpdt_97*[1 1],'-.k','LineWidth',lineWidth)
p(1) = plot(out.t([1 end]),1e-3*dpdt_97*[-1 -1],'-.k','LineWidth',lineWidth);
p(1).HandleVisibility='off';

dist_dpdt = statsTimeVar_cdf(out.t,abs(out.dydt(:,par.iy.p_ro)));
dpdt_99 = dist_dpdt.xi(find(dist_dpdt.f > 0.99,1,'first'));
plot(out.t([1 end]),1e-3*dpdt_99*[1 1],':k','LineWidth',lineWidth*2)
p(2) = plot(out.t([1 end]),1e-3*dpdt_99*[-1 -1],':k','LineWidth',lineWidth*2);
p(2).HandleVisibility='off';

xlabel('time (s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('rate (kPa/s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

% yLim = ylim;
% ylim([0 yLim(2)])

leg = legend('$dp_{ro}/dt$','target limit','$\pm P_{97} |dp_{ro}/dt|$', ...
             '$\pm P_{99} |dp_{ro}/dt|$', ...
             'Interpreter','latex');
leg.FontSize = axFontSize;
leg.FontName = 'Times';
leg.Orientation = 'horizontal';
pos = leg.Position;
leg.Position = [0.5-pos(3)/2, 0, pos(3), pos(4)];
leg.Box = 'off';

title(titleString,...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

ax = gca;
ax.FontName = 'times';
ax.FontSize = axFontSize;
pos = ax.Position;
adj = 0.125;
ax.Position = [pos(1), pos(2)+adj, pos(3), pos(4)-adj];