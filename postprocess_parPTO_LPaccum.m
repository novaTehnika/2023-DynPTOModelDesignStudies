%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_LPaccum")
        load(files(j).name)

        % postprocess simulation data
        if sum(imag(out.p_lout+out.p_a+out.p_b)) == 0

             % Mean pressure at WEC-driven pump inlet
            p_loutMean(iVar) = mean(out.p_lout);
             % Variation in pressure at WEC-driven pump inlet
            p_loutMax(iVar) = max(out.p_lout);
            p_loutMin(iVar) = min(out.p_lout);
            p_loutVar(iVar) = var(out.p_lout);
            p_loutStd(iVar) = std(out.p_lout);
             % Minimum pressure in WEC-driven pump chambers
            p_wpMin(iVar) = min(min(out.p_a),min(out.p_b));
    
             % Electric power consumption of charge pump
            P_cElec(iVar) = mean(out.power.P_cElec);
            P_cElec_norm(iVar) = P_cElec(iVar)/mean(out.power.P_WEC);
             % Power losses from charge pump
            P_cLoss(iVar) = mean(out.power.P_cLoss);
            L_c(iVar) = P_cLoss(iVar)/mean(out.power.P_WEC);

             % power loss from pipeline
            P_LPPL(iVar) = out.power.P_LPPL;
            L_LPPL(iVar) = P_LPPL(iVar)/mean(out.power.P_WEC);

        else

            p_loutMean(iVar) = nan;
             % Variation in pressure at WEC-driven pump inlet
            p_loutMax(iVar) = nan;
            p_loutMin(iVar) = nan;
            p_loutVar(iVar) = nan;
            p_loutStd(iVar) = nan;
             % Minimum pressure in WEC-driven pump chambers
            p_wpMin(iVar) = nan;
    
             % Electric power consumption of charge pump
            P_cElec(iVar) = nan;
            P_cElec_norm(iVar) = nan;
             % Power losses from charge pump
            P_cLoss(iVar) = nan;
            L_c(iVar) = nan;

             % power loss from pipeline
            P_LPPL(iVar) = nan;
            L_LPPL(iVar) = nan;

        end

    end

end

clearvars out

%% Find optimal distribution of accumulator volume for each total 
% accumulator volume and pipeline diameter
 % find individuals meeting cavitation constraints
p_cavLimit = 0.5e4; % [Pa] prescribed limit on pressure in WEC-driven pump
[~,meetsConstraints] = find(p_wpMin >= p_cavLimit);
 % find non-dominated individuals from set meeting dpdt criterion
non_dominated = paretoFront2D(Vc_l_mesh(meetsConstraints),'min', ...
               P_LPPL(meetsOCnstraints) + P_cElec(meetsConstraints),'min');
[~, ii_sort] = sort(Vc_l_mesh(meetsConstraints(non_dominated)));
iiPareto = meetsConstraints(non_dominated(ii_sort));
Vc_l_opt = Vc_l_mesh(iiPareto);
X_opt = w_c_mesh(iiPareto);
d_LPPL = d_LPPL_mesh(iiPareto);
clearvars meetsConstraints non_dominated ii_sort


%% Plot Pereto optimal results

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
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 10;
fontSize = 9;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

titleString = ['Performance of Low-Pressure Circuit Branch',newline,...
                'as a Function of Installed Low-Pressure Accumulator Volume',newline,...
                'Sea State ',num2str(SS)];
sgtitle(titleString,...
'Interpreter','latex','FontSize',fontSize+2,'fontname','Times')

n_plots = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure WEC-driven pump inlet
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,1e-5*p_loutMean(iiPareto),'k','Marker','x');
ip = ip+1;

p(ip,iax) = plot(Vc_l_opt,1e-5*p_loutMax(iiPareto),'r','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,1e-5*p_loutMin(iiPareto),'r','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

p(ip,iax) = plot(Vc_l_opt,1e-5*(p_loutMean(iiPareto)+p_loutStd(iiPareto)),':k','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,1e-5*(p_loutMean(iiPareto)-p_loutStd(iiPareto)),':k','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('pressure (bar)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Pressure at WEC-Driven Pump Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend('mean','max & min','mean +- 1 stdDev');
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xLim = xlim;
xlim([0 xLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure in WEC-driven pump
iax = 2;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt([1 end]),1e-5*p_cavLimit*[1 1],'r');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,1e-5*p_wpMin(iiPareto),'k','Marker','x');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('pressure (bar)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Minimum Pressure in WEC-Driven Pump'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend('limit');
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power Consumption of Charge Pump
iax = 3;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,1e-3*P_cElec(iiPareto),'k','Marker','x');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power (kW)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Mean Electric Power Consumption of Charge Pump'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power loss
iax = 4;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,100*L_c(iiPareto),'r','Marker','^');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,100*L_LPPL(iiPareto),'b','Marker','o');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt, ...
                100*(L_c(iiPareto)+L_LPPL(iiPareto)),'k','Marker','x');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Mean Power Loss of Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend('charge pump','pipeline','combined');
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LP pipeline diameter
iax = 5;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,100*d_LPPL_opt,'k','Marker','x');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('diameter (cm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Internal Diameter of the Low-Pressure Pipeline'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%% Plot all results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 10;
fontSize = 9;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

titleString = ['Performance of Low-Pressure Circuit Branch',newline,...
                'as a Function of Installed Low-Pressure Accumulator Volume',newline,...
                'Sea State ',num2str(SS)];
sgtitle(titleString,...
'Interpreter','latex','FontSize',fontSize+2,'fontname','Times')

n_plots = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure WEC-driven pump inlet
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = scatter(Vc_l_mesh,1e-5*p_loutMean,'xk');
ip = ip+1;

p(ip,iax) = scatter(Vc_l_mesh,1e-5*p_loutMax,'xr');
ip = ip+1;
p(ip,iax) = scatter(Vc_l_mesh,1e-5*p_loutMin,'xr');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

p(ip,iax) = scatter(Vc_l_mesh,1e-5*(p_loutMean+p_loutStd),'.k');
ip = ip+1;
p(ip,iax) = scatter(Vc_l_mesh,1e-5*(p_loutMean-p_loutStd),'.k');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('pressure (bar)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Pressure at RO Module Feed Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend('mean','max & min','mean +- 1 stdDev');
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xLim = xlim;
xlim([0 xLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure in WEC-driven pump
iax = 2;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(Vc_l([1 end]),1e-5*p_cavLimit*[1 1],'r');
ip = ip+1;
p(ip,iax) = scatter(Vc_l_mesh,1e-5*p_wpMin,'xk');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('pressure (bar)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Minimum Pressure in WEC-Driven Pump'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend('limit');
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power Consumption of CHarge Pump
iax = 3;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = scatter(Vc_l_mesh,1e-3*P_cElec,'xk');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power (kW)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Mean Electric Power Consumption of Charge Pump'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power loss
iax = 4;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = scatter(Vc_l_mesh,100*L_c,'xk');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Mean Power Loss of Charge Pump Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Charge Pump Shaft Speed
iax = 5;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = scatter(Vc_l_mesh,60/2/pi*w_c_mesh,'xk');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('speed (rpm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Mean Power Loss of Charge Pump Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])