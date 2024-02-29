%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_motorGen")
        load(files(j).name)

        % postprocess simulation data
         % Mean pressure at RO module
        p_roMean(iVar) = mean(out.p_ro);
         % Variation in pressure at RO module
        p_roMax(iVar) = max(out.p_ro);
        p_roMin(iVar) = min(out.p_ro);
        p_roVar(iVar) = var(out.p_ro);
        p_roStd(iVar) = std(out.p_ro);
         % Rate of change in pressure at RO module
          % max
        dpdt_max(iVar) = max(out.dydt(:,out.par.iy.p_ro));
          % 97th percentile
        dist_dpdt = statsTimeVar_cdf(out.t,abs(out.dydt(:,out.par.iy.p_ro)));
        dpdt_97(iVar) = dist_dpdt.xi(find(dist_dpdt.f > 0.97,1,'first'));
        clearvars dist_dpdt

         % Electric power generation
        P_gen(iVar) = mean(out.power.P_gen);
         % Power losses from motor and generator
        P_genLoss(iVar) = mean(out.power.P_genLoss);
        P_pmLoss(iVar) = mean(out.power.P_pmLoss);
        L_pm(iVar) = P_pmLoss(iVar)/mean(out.power.P_WEC);
        L_gen(iVar) = P_genLoss(iVar)/mean(out.power.P_WEC);

         % Permeate production
        q_perm(iVar) = mean(out.q_perm);

         % Power losses from pressure relief valves
        P_prv(iVar) = mean(out.power.P_hinPRV + out.power.P_roPRV);
        L_prv(iVar) = P_prv(iVar)/mean(out.power.P_WEC);


    end

end

clearvars out

%% Plot average power loss from the ripple control valve as a function of total accumulator volume for distribution (color) and valve coefficient (line type)
DpmUnitConvert = 1/(1e-6/(2*pi)); % [m^3/rad -> cc/rev]

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

supTitleFontSize = 10;
subTitleFontSize = 9;
axFontSize = 8;
bottomEdge = 1;
leftEdge = 3;
width = 5.75; % one column: 3+9/16, two column: 7.5
height = 7;
lineWidth = 0.5;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

switch par.rvConfig.active
    case 0
        rvStr = 'Passive Ripple Control';
    case 1
        rvStr = 'Active Ripple Control';
end
titleString = ['Hydraulic Motor Sizing Study with ',rvStr,':',newline,...
                'Sea State ',num2str(SS)];
sgtitle(titleString,...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

n_plots = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure
iax = 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,1e-6*p_roMean,'k','Marker','x');
ip = ip+1;

p(ip,iax) = plot(Dpm*DpmUnitConvert,1e-6*p_roMax,'r','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,1e-6*p_roMin,'r','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

p(ip,iax) = plot(Dpm*DpmUnitConvert,1e-6*(p_roMean+p_roStd),':k','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,1e-6*(p_roMean-p_roStd),':k','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

grid on

xlabel('displacement (cc/rev) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('pressure (MPa)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Pressure at RO Module Feed Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

axL = gca;
axL.FontName = 'Liberation Serif';
axL.FontSize = axFontSize;

leg = legend('mean','max. and min.','mean +- 1 stdDev');
leg.FontSize = axFontSize;
leg.Interpreter = 'latex';
leg.FontName = 'Times';
set(leg, 'Location', 'best')


xLim = xlim;
xlim([0 xLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rate of change in pressure
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Dpm([1 end])*DpmUnitConvert,1e-3*par.control.dpdt_ROmax*[1 1],'r','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,1e-3*dpdt_max,'k','Marker','x');
ip = ip+1;

p(ip,iax) = plot(Dpm*DpmUnitConvert,1e-3*dpdt_97,'--k','Marker','x');
ip = ip+1;

grid on

xlabel('displacement (cc/rev) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('rate of change (kPa/s)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Rate of Change in Pressure at RO Module Feed Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

axL = gca;
axL.FontName = 'Liberation Serif';
axL.FontSize = axFontSize;

leg = legend('limit','max','P97');
leg.FontSize = axFontSize;
leg.Interpreter = 'latex';
leg.FontName = 'Times';
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power generation
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,1e-3*P_gen,'k','Marker','x');
ip = ip+1;

grid on

xlabel('displacement (cc/rev) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('power (kW)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Mean Electric Power Generation'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

axL = gca;
axL.FontName = 'Liberation Serif';
axL.FontSize = axFontSize;

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power loss
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,L_pm,'k','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,L_gen,'--k','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,L_prv,'-.k','Marker','x');
ip = ip+1;

grid on

xlabel('displacement (cc/rev) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('power loss', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Mean Power Loss Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

axL = gca;
axL.FontName = 'Liberation Serif';
axL.FontSize = axFontSize;

leg = legend('hyd. motor','generator','pressure relief valves');
leg.FontSize = axFontSize;
leg.Interpreter = 'latex';
leg.FontName = 'Times';
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Permeate production
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Dpm*DpmUnitConvert,24*3600*q_perm,'k','Marker','x');
ip = ip+1;

grid on

xlabel('displacement (cc/rev) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('production rate (m\textsuperscript{3}/day)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Mean Rate of Permeate Production'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

axL = gca;
axL.FontName = 'Liberation Serif';
axL.FontSize = axFontSize;

xLim = xlim;
xlim([0 xLim(2)])
