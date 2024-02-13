%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_timeSpanConvergence")
        load(files(j).name,'-regexp','^(?!out)\w')

        % Calculate metrics
        studyData.Ebal_error(iVar,SS) = Ebal_error;
        studyData.Vbal_error(iVar,SS) = Vbal_error;
        
        % for electrical power balance
        studyData.deltaE_battery(iVar,SS) = deltaE_battery;
        
        % for LPaccum
        studyData.p_loutMean(iVar,SS) = p_loutMean;
         % Variation in pressure at WEC-driven pump inlet
        studyData.p_loutMax(iVar,SS) = p_loutMax;
        studyData.p_loutMin(iVar,SS) = p_loutMin;
        studyData.p_loutVar(iVar,SS) = p_loutVar;
        studyData.p_loutStd(iVar,SS) = p_loutStd;
         % Minimum pressure in WEC-driven pump chambers
        studyData.p_wpMin(iVar,SS) = p_wpMin;
        
         % Electric power consumption of charge pump
        studyData.P_cElec(iVar,SS) = P_cElec;
        studyData.P_cElec_norm(iVar,SS) = P_cElec_norm;
         % Power losses from charge pump
        studyData.P_cLoss(iVar,SS) = P_cLoss;
        studyData.L_c(iVar,SS) = L_c;
        
         % power loss from pipeline
        studyData.P_LPPL(iVar,SS) = mean(P_LPPL);
        studyData.L_LPPL(iVar,SS) = mean(L_LPPL);
        
        % for accum_woRV and accum_wRV
        studyData.q_permMean(iVar,SS) = q_permMean;
        studyData.PP_WEC(iVar,SS) = PP_WEC;
        studyData.PP_wp(iVar,SS) = PP_wp;
        studyData.PP_rv(iVar,SS) = PP_rv;
        studyData.PP_hinPRV(iVar,SS) = PP_hinPRV;
        studyData.PP_roPRV(iVar,SS) = PP_roPRV;
        studyData.dpdt_max(iVar,SS) = dpdt_max;
        studyData.dpdt_97(iVar,SS) = dpdt_97;

         % power loss from pipeline
        studyData.P_HPPL(iVar,SS) = mean(P_HPPL);
        studyData.L_HPPL(iVar,SS) = mean(L_HPPL);

    end

end

%% Find indices for missing data files

% ...

%% Plotting
SS = 3;
err_endNorm = @(x) x./x(end); 

black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];

linestyles = {'-', '--', ':', '-.'};

%% Plot energy and volume balance error as a function of timestep
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

xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Normalized Energy and Volume Balance: Sea State',num2str(SS)], ...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')
yyaxis left
plot(tspan,abs(studyData.Ebal_error(:,SS)));
hold on
ylabel('energy balance error, normalized', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')


yyaxis right
hold on
plot(tspan,abs(studyData.Vbal_error(:,SS)));
ylabel('volume balance error, normalized', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

%% 
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
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

% for electrical power balance
plot(tspan,studyData.deltaE_battery(:,SS));



%% for LPaccum
 % Variation in pressure at WEC-driven pump inlet
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 8;
fontSize = 9;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 5;
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.p_loutMean(:,SS)));
plot(tspan,err_endNorm(studyData.p_loutMax(:,SS)));
plot(tspan,err_endNorm(studyData.p_loutMin(:,SS)));

legend('mean(p_{lout})','max(p_{lout})','min(p_{lout})')

iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.p_loutVar(:,SS)));
plot(tspan,err_endNorm(studyData.p_loutStd(:,SS)));

legend('var(p_{lout})','(stddev(p_{lout})')

 % Minimum pressure in WEC-driven pump chambers
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.p_wpMin(:,SS)));

legend('min(p_{wp})')

 % Electric power consumption of charge pump &
 % Power losses from charge pump
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.P_cElec(:,SS)));
plot(tspan,err_endNorm(studyData.P_cElec_norm(:,SS)));
plot(tspan,err_endNorm(studyData.P_cLoss(:,SS)));
plot(tspan,err_endNorm(studyData.L_c(:,SS)));

legend('P_cElec','P_cElec_norm','P_cLoss','L_c')

 % power loss from pipeline
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.P_LPPL(:,SS)));
plot(tspan,err_endNorm(studyData.L_LPPL(:,SS)));

legend('P_LPPL','L_LPPL')

%% permeate
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
width = 5.75; % one column: 3+9/16, two column: 7.5
height = 2;
supTitleFontSize = 9;
subTitleFontSize = 9;
axFontSize = 8;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

semilogy(tspan,err_endNorm(studyData.q_permMean(:,SS)),'color',black,'LineWidth',lineWidth)

grid on

ylabel('error', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
title(['Permeate Production Convergence Study: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

% leg = legend('energy','volume', ...
%              'Interpreter','latex');
% leg.FontSize = axFontSize;
% leg.FontName = 'Times';
% leg.Orientation = 'horizontal';
% pos = leg.Position;
% leg.Position = [0.5-pos(3)/2, 0, pos(3), pos(4)];
% leg.Box = 'off';

ax = gca;
ax.FontName = 'Times';
ax.FontSize = axFontSize;

%% power
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
width = 5.75; % one column: 3+9/16, two column: 7.5
height = 2;
supTitleFontSize = 9;
subTitleFontSize = 9;
axFontSize = 8;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];


semilogy(tspan,err_endNorm(studyData.PP_WEC(:,SS)),'color',black,'LineWidth',lineWidth)

grid on

ylabel('error', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
title(['WEC Power Capture Convergence Study: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

% leg = legend('energy','volume', ...
%              'Interpreter','latex');
% leg.FontSize = axFontSize;
% leg.FontName = 'Times';
% leg.Orientation = 'horizontal';
% pos = leg.Position;
% leg.Position = [0.5-pos(3)/2, 0, pos(3), pos(4)];
% leg.Box = 'off';

ax = gca;
ax.FontName = 'Times';
ax.FontSize = axFontSize;

%% pressure dpdt
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
width = 5.75; % one column: 3+9/16, two column: 7.5
height = 2;
supTitleFontSize = 9;
subTitleFontSize = 9;
axFontSize = 8;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

plot(tspan,err_endNorm(studyData.dpdt_max(:,SS)));
plot(tspan,err_endNorm(studyData.dpdt_97(:,SS)));

semilogy(tspan,err_endNorm(studyData.dpdt_max(:,SS)),'color',black,'LineWidth',lineWidth)
hold on
semilogy(tspan,err_endNorm(studyData.dpdt_97(:,SS)),'color',maroon,'LineWidth',lineWidth)
grid on

ylabel('error', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
title(['Rate of Change in Feed Pressure Convergence Study: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

leg = legend('$|dp_f/dt|$','$P_{97}|dp_f/dt|$', ...
             'Interpreter','latex');
leg.FontSize = axFontSize;
leg.FontName = 'Times';
% leg.Orientation = 'horizontal';
% pos = leg.Position;
% leg.Position = [0.5-pos(3)/2, 0, pos(3), pos(4)];
% leg.Box = 'off';

ax = gca;
ax.FontName = 'Times';
ax.FontSize = axFontSize;

%% for accum_woRV and accum_wRV
bottomEdge = 1;
leftEdge = 3;
width = 5.75; % one column: 3+9/16, two column: 7.5
height = 8;
fontSize = 9;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 4;
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.q_permMean(:,SS)));

legend('q\_permMean')

iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.PP_WEC(:,SS)));
plot(tspan,err_endNorm(studyData.PP_wp(:,SS)));
plot(tspan,err_endNorm(studyData.PP_rv(:,SS)));
plot(tspan,err_endNorm(studyData.PP_hinPRV(:,SS)));
plot(tspan,err_endNorm(studyData.PP_roPRV(:,SS)));

legend('PP_WEC','PP_wp','PP_rv','PP_hinPRV','PP_roPRV')

iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.dpdt_max(:,SS)));
plot(tspan,err_endNorm(studyData.dpdt_97(:,SS)));

legend('dpdt_max','dpdt_97')

 % power loss from pipeline
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;
xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

hold on

plot(tspan,err_endNorm(studyData.P_HPPL(:,SS)));
plot(tspan,err_endNorm(studyData.L_HPPL(:,SS)));

legend('P_HPPL','L_HPPL')

xlabel('time span (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

