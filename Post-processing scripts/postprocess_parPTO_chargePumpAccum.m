%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_chargePumpAccum")
        load(files(j).name)

        % save post-processed data

         % Mean pressure at WEC-driven pump inlet
        studyData.p_loutMean(iVar) = p_loutMean;
         % Variation in pressure at WEC-driven pump inlet
        studyData.p_loutMax(iVar) = p_loutMax;
        studyData.p_loutMin(iVar) = p_loutMin;
        studyData.p_loutVar(iVar) = p_loutVar;
        studyData.p_loutStd(iVar) = p_loutStd;
         % Minimum pressure in WEC-driven pump chambers
        studyData.p_wpMin(iVar) = p_wpMin;

         % Electric power consumption of charge pump
        studyData.P_cElec(iVar) = P_cElec;
        studyData.L_cElec(iVar) = L_cElec;
         % Power losses from charge pump (inlcudes electric motor loss)
        studyData.P_cLoss(iVar) = P_cLoss;
        studyData.L_c(iVar) = L_c;

    end

end

clearvars out

%%



%% Find indices for missing data files

files = dir;
nfiles = size(files,1);
notDone = 1:nVar;
Done = [];


for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_chargePumpAccum")
        load(files(j).name,'SS','iVar')
        [r,c,val] = find(notDone==iVar);
        notDone = [notDone(1:c-1), notDone(c+1:end)];
        Done = [Done, iVar];

    end

end


    try
        doneArrayStr = num2str(Done(1));
        for j = 2:length(Done)
            switch 1
                case 1
                    doneArrayStr = append(doneArrayStr,[',',num2str(Done(j))]);
                case 2
                    doneArrayStr = append(doneArrayStr,[',',num2str(Done(j),['%0',floor(log10(nVar)),'.f'])]);
            end
        end
    catch
        % just move o
    end

    try
        jobArrayStr = num2str(notDone(1));
        for j = 2:length(notDone)
            jobArrayStr = append(jobArrayStr,[',',num2str(notDone(j))]);
        end

        if 1
            for j = 1:length(notDone)
                iVar = notDone(j);
                studyData.p_loutMean(iVar) = nan;
                 % Variation in pressure at WEC-driven pump inlet
                studyData.p_loutMax(iVar) = nan;
                studyData.p_loutMin(iVar) = nan;
                studyData.p_loutVar(iVar) = nan;
                studyData.p_loutStd(iVar) = nan;
                 % Minimum pressure in WEC-driven pump chambers
                studyData.p_wpMin(iVar) = nan;
        
                 % Electric power consumption of charge pump
                studyData.P_cElec(iVar) = nan;
                studyData.L_cElec(iVar) = nan;
                 % Power losses from charge pump
                studyData.P_cLoss(iVar) = nan;
                studyData.L_c(iVar) = nan;

            end
        end

    catch
        % just move on
    end

%% Find optimal charge pump speed for each total accumulator volume
 % find individuals meeting cavitation constraints
p_cavLimit = 0.5e4; % [Pa] prescribed limit on pressure in WEC-driven pump
[~,meetsConstraints] = find(studyData.p_wpMin >= p_cavLimit);
 % find non-dominated individuals from set meeting dpdt criterion
non_dominated = paretoFront2D(Vc_l_mesh(meetsConstraints),'min', ...
                              studyData.L_c(meetsConstraints),'min');
[~, ii_sort] = sort(Vc_l_mesh(meetsConstraints(non_dominated)));
iiPareto = meetsConstraints(non_dominated(ii_sort));
Vc_l_opt = Vc_l_mesh(iiPareto);
w_c_opt = w_c_mesh(iiPareto);
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

titleString = ['Performance of Low-Pressure Circuit Branch',newline,...
                'as a Function of Installed Low-Pressure Accumulator Volume',newline,...
                'Sea State ',num2str(SS)];
sgtitle(titleString,...
'Interpreter','latex','FontSize',supTitleFontSize,'fontname','Times')

n_plots = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure WEC-driven pump inlet
iax = 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,1e-5*studyData.p_loutMean(iiPareto),'k','Marker','x');
ip = ip+1;

p(ip,iax) = plot(Vc_l_opt,1e-5*studyData.p_loutMax(iiPareto),'r','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,1e-5*studyData.p_loutMin(iiPareto),'r','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

p(ip,iax) = plot(Vc_l_opt,1e-5*(studyData.p_loutMean(iiPareto)+studyData.p_loutStd(iiPareto)),':k','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,1e-5*(studyData.p_loutMean(iiPareto)-studyData.p_loutStd(iiPareto)),':k','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('pressure (bar)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Pressure at RO Module Feed Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

leg = legend('mean','max {\&} min','mean +- 1 stdDev');
leg.FontSize = axFontSize;
leg.Interpreter = 'latex';
leg.FontName = 'Times';
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 xLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure in WEC-driven pump
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt([1 end]),1e-5*p_cavLimit*[1 1],'r');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,1e-5*studyData.p_wpMin(iiPareto),'k','Marker','x');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('pressure (bar)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Minimum Pressure in WEC-Driven Pump'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

axL = gca;
axL.FontName = 'Liberation Serif';
axL.FontSize = axFontSize;

leg = legend('limit');
leg.FontSize = axFontSize;
leg.Interpreter = 'latex';
leg.FontName = 'Times';
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power Consumption of Charge Pump
% iax = iax + 1;
% ax(iax) = subplot(n_plots,1,iax);
% 
% hold on
% 
% ip = 1;
% p(ip,iax) = plot(Vc_l_opt,1e-3*studyData.P_cElec(iiPareto),'k','Marker','x');
% ip = ip+1;
% 
% xlabel('volume (1000L) ', ...
% 'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
% ylabel('power (kW)', ...
% 'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
% 
% titleString = ['Mean Electric Power Consumption of Charge Pump'];
% title(titleString,...
% 'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')
% 
% axL = gca;
% axL.FontName = 'Liberation Serif';
% axL.FontSize = axFontSize;
%
% xLim = xlim;
% xlim([0 xLim(2)])
% yLim = ylim;
% ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power loss
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,100*studyData.L_c(iiPareto),'k','Marker','x');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Mean Power Loss of Charge Pump Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

axL = gca;
axL.FontName = 'Liberation Serif';
axL.FontSize = axFontSize;

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Charge Pump Shaft Speed
iax = iax + 1;
ax(iax) = subplot(n_plots,1,iax);

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,60/2/pi*w_c_opt,'k','Marker','x');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')
ylabel('speed (rpm)', ...
'Interpreter','latex','FontSize',axFontSize,'fontname','Times')

titleString = ['Mean Power Loss of Charge Pump Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',subTitleFontSize,'fontname','Times')

axL = gca;
axL.FontName = 'Liberation Serif';
axL.FontSize = axFontSize;

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
p(ip,iax) = scatter(Vc_l_mesh,1e-5*studyData.p_loutMean,'xk');
ip = ip+1;

p(ip,iax) = scatter(Vc_l_mesh,1e-5*studyData.p_loutMax,'xr');
ip = ip+1;
p(ip,iax) = scatter(Vc_l_mesh,1e-5*studyData.p_loutMin,'xr');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

p(ip,iax) = scatter(Vc_l_mesh,1e-5*(studyData.p_loutMean+studyData.p_loutStd),'.k');
ip = ip+1;
p(ip,iax) = scatter(Vc_l_mesh,1e-5*(studyData.p_loutMean-studyData.p_loutStd),'.k');
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
p(ip,iax) = scatter(Vc_l_mesh,1e-5*studyData.p_wpMin,'xk');
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
p(ip,iax) = scatter(Vc_l_mesh,1e-3*studyData.P_cElec,'xk');
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
p(ip,iax) = scatter(Vc_l_mesh,100*studyData.L_c,'xk');
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