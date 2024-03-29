%% Collect data from data files
files = dir;
nfiles = size(files,1);
dataFileName = "data_parPTO_LPaccum";
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,dataFileName)
        load(files(j).name)

        for iw_c = 1:nVar4

            % study variables
            studyVar(SS).Vc_l(iVar,iw_c) = Vc_l_mesh(iVar);
            studyVar(SS).X(iVar,iw_c) = X_mesh(iVar);
            studyVar(SS).d_line(iVar,iw_c) = d_LPPL_mesh(iVar);
            studyVar(SS).w_c(iVar,iw_c) = w_c(iw_c);

            % metrics
             % Mean pressure at WEC-driven pump inlet
            studyData(SS).p_loutMean(iVar,iw_c) = p_loutMean(iw_c);
             % Variation in pressure at WEC-driven pump inlet
            studyData(SS).p_loutMax(iVar,iw_c) = p_loutMax(iw_c);
            studyData(SS).p_loutMin(iVar,iw_c) = p_loutMin(iw_c);
            studyData(SS).p_loutVar(iVar,iw_c) = p_loutVar(iw_c);
            studyData(SS).p_loutStd(iVar,iw_c) = p_loutStd(iw_c);
             % Minimum pressure in WEC-driven pump chambers
            studyData(SS).p_wpMin(iVar,iw_c) = p_wpMin(iw_c);
            studyData(SS).p_staticCVmin(iVar,iw_c) = p_staticCVmin(iw_c);
    
             % Electric power consumption of charge pump
            studyData(SS).P_cElec(iVar,iw_c) = P_cElec(iw_c);
            studyData(SS).P_cElec_norm(iVar,iw_c) = P_cElec_norm(iw_c);
             % Power losses from charge pump
            studyData(SS).P_cLoss(iVar,iw_c) = P_cLoss(iw_c);
            studyData(SS).L_c(iVar,iw_c) = L_c(iw_c);

             % power loss from pipeline
            studyData(SS).P_LPPL(iVar,iw_c) = P_LPPL(iw_c);
            studyData(SS).L_LPPL(iVar,iw_c) = L_LPPL(iw_c);

        end

    end

end

%% Find indices for missing data files

files = dir;
nfiles = size(files,1);
for SS = 1:6
    notDone(SS).list = 1:nVar;
    Done(SS).list = [];
end

for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,dataFileName)
        load(files(j).name,'SS','iVar')
        [r,c,val] = find(notDone(SS).list==iVar);
        notDone(SS).list = [notDone(SS).list(1:c-1), notDone(SS).list(c+1:end)];
        Done(SS).list = [Done(SS).list, iVar];

    end

end

for SS = 1:6
    try
        doneArrayStr(SS).string = num2str(Done(SS).list(1));
        for j = 2:length(Done(SS).list)
            switch 1
                case 1
                    doneArrayStr(SS).string = append(doneArrayStr(SS).string,[',',num2str(Done(SS).list(j))]);
                case 2
                    doneArrayStr(SS).string = append(doneArrayStr(SS).string,[',',num2str(Done(SS).list(j),['%0',floor(log10(nVar)),'.f'])]);
            end
        end
    catch
        % just move on
    end

    try
        jobArrayStr(SS).string = num2str(notDone(SS).list(1));
        for j = 2:length(notDone(SS).list)
            jobArrayStr(SS).string = append(jobArrayStr(SS).string,[',',num2str(notDone(SS).list(j))]);
        end

        if 1
            nanArray = nan*ones(nVar4,1);
            for j = 1:length(notDone(SS).list)
                iVar = notDone(SS).list(j);
                % study variables
                studyVar(SS).Vc_l(iVar,:) = ones(1,nVar4)*Vc_l_mesh(iVar);
                studyVar(SS).X(iVar,:) = ones(1,nVar4)*X_mesh(iVar);
                studyVar(SS).d_line(iVar,:) = ones(1,nVar4)*d_LPPL_mesh(iVar);
                studyVar(SS).w_c(iVar,:) = w_c;

                 % Mean pressure at WEC-driven pump inlet
                studyData(SS).p_loutMean(iVar,:) = nanArray;
                 % Variation in pressure at WEC-driven pump inlet
                studyData(SS).p_loutMax(iVar,:) = nanArray;
                studyData(SS).p_loutMin(iVar,:) = nanArray;
                studyData(SS).p_loutVar(iVar,:) = nanArray;
                studyData(SS).p_loutStd(iVar,:) = nanArray;
                 % Minimum pressure in WEC-driven pump chambers
                studyData(SS).p_wpMin(iVar,:) = nanArray;
                studyData(SS).p_staticCVmin(iVar,:) = nanArray;

                 % Electric power consumption of charge pump
                studyData(SS).P_cElec(iVar,:) = nanArray;
                studyData(SS).P_cElec_norm(iVar,:) = nanArray;
                 % Power losses from charge pump
                studyData(SS).P_cLoss(iVar,:) = nanArray;
                studyData(SS).L_c(iVar,:) = nanArray;

                 % power loss from pipeline
                studyData(SS).P_LPPL(iVar,:) = nanArray;
                studyData(SS).L_LPPL(iVar,:) = nanArray;

            end
        end

    catch
        % just move on
    end
end

%% Find optimal distribution of accumulator volume for each total 
SS = 2;

% accumulator volume and pipeline diameter
 % find individuals meeting cavitation constraints
p_cavLimit = 10e3; % [Pa] prescribed limit on pressure in WEC-driven pump
id_line = 10;
d_line_nom = d_LPPL(id_line);
% meetsConstraints = find(studyData(SS).p_wpMin(:) >= p_cavLimit ...
%                         & studyVar(SS).d_line(:) == d_line_nom);
meetsConstraints = find(studyData(SS).p_staticCVmin(:) >= p_cavLimit ...
                        & studyVar(SS).d_line(:) == d_line_nom);

 % find non-dominated individuals from set meeting cavitation constraints
non_dominated = paretoFront2D(studyVar(SS).Vc_l(meetsConstraints),'min', ...
               studyData(SS).P_LPPL(meetsConstraints) ...
               + studyData(SS).P_cElec(meetsConstraints),'min');
[~, ii_sort] = sort(studyVar(SS).Vc_l(meetsConstraints(non_dominated)));
iiPareto = meetsConstraints(non_dominated(ii_sort));

Vc_l_opt = studyVar(SS).Vc_l(iiPareto);
Vc_l_meetsCon = studyVar(SS).Vc_l(meetsConstraints);

X_opt = studyVar(SS).X(iiPareto);
X_meetsCon = studyVar(SS).X(meetsConstraints);

d_LPPL_opt = studyVar(SS).d_line(iiPareto);

w_c_opt = studyVar(SS).w_c(iiPareto);
w_c_meetsCon = studyVar(SS).w_c(meetsConstraints);

% clearvars meetsConstraints non_dominated ii_sort

%% Plot Pareto optimal results

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
                'Sea State ',num2str(SS),newline,...
                'Pipeline Diameter of ',num2str(100*d_line_nom),' cm'];
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
p(ip,iax) = plot(Vc_l_opt,1e-5*studyData(SS).p_loutMean(iiPareto),'k','Marker','x');
ip = ip+1;

p(ip,iax) = plot(Vc_l_opt,1e-5*studyData(SS).p_loutMax(iiPareto),'r','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,1e-5*studyData(SS).p_loutMin(iiPareto),'r','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

p(ip,iax) = plot(Vc_l_opt,1e-5*(studyData(SS).p_loutMean(iiPareto)+studyData(SS).p_loutStd(iiPareto)),':k','Marker','x');
ip = ip+1;
p(ip,iax) = plot(Vc_l_opt,1e-5*(studyData(SS).p_loutMean(iiPareto)-studyData(SS).p_loutStd(iiPareto)),':k','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

grid on

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
p(ip,iax) = plot(Vc_l_opt,1e-5*studyData(SS).p_wpMin(iiPareto),'k','Marker','x');
ip = ip+1;

grid on

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power Consumption of Charge Pump
% iax = 3;
% ax(iax) = subplot(n_plots,1,iax);
% ax(iax).FontName = 'times';
% ax(iax).FontSize = fontSize-1;
% 
% hold on
% 
% ip = 1;
% p(ip,iax) = plot(Vc_l_opt,1e-3*studyData(SS).P_cElec(iiPareto),'k','Marker','x');
% ip = ip+1;
% 
% grid on
% 
% xlabel('volume (1000L) ', ...
% 'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
% ylabel('power (kW)', ...
% 'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
% 
% titleString = ['Mean Electric Power Consumption of Charge Pump'];
% title(titleString,...
% 'Interpreter','latex','FontSize',fontSize,'fontname','Times')
% 
% xLim = xlim;
% xlim([0 xLim(2)])
% yLim = ylim;
% ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power loss
iax = 3;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;



ip = 1;
p(ip,iax) = semilogy(Vc_l_opt,100*studyData(SS).L_c(iiPareto),'r','Marker','^');
hold on

ip = ip+1;
p(ip,iax) = semilogy(Vc_l_opt,100*studyData(SS).L_LPPL(iiPareto),'b','Marker','o');
ip = ip+1;
p(ip,iax) = semilogy(Vc_l_opt, ...
                100*(studyData(SS).L_c(iiPareto)+studyData(SS).L_LPPL(iiPareto)),'k','Marker','x');
ip = ip+1;

grid on

yLim = ylim;
ylim([yLim(1) 100])

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% charge pump shaft speed
iax = 4;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,w_c_opt/2/pi*60,'k','Marker','x');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('speed (rpm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Charge Pump Shaft Speed'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% distibution of volume
iax = 5;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(Vc_l_opt,X_opt,'k','Marker','x');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('portion', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Portion of Volume at WEC-driven Pump Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%% Plot all results meeting constraints

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
                'Sea State ',num2str(SS),newline,...
                'Pipeline Diameter of ',num2str(100*d_line_nom),' cm'];
sgtitle(titleString,...
'Interpreter','latex','FontSize',fontSize+2,'fontname','Times')

n_plots = 6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pressure WEC-driven pump inlet
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = scatter(Vc_l_meetsCon,1e-5*studyData(SS).p_loutMin(meetsConstraints),'r','Marker','x');
p(ip,iax).HandleVisibility='off';
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('pressure (bar)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Min. Pressure at WEC-Driven Pump Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

% set(leg, 'Position', rect)
% set(leg, 'Location', 'best')
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
p(ip,iax) = plot([min(Vc_l_meetsCon) max(Vc_l_meetsCon)],1e-5*p_cavLimit*[1 1],'r');
ip = ip+1;
% p(ip,iax) = scatter(Vc_l_meetsCon,1e-5*studyData(SS).p_wpMin(meetsConstraints),'k','Marker','x');
p(ip,iax) = scatter(Vc_l_meetsCon,1e-5*studyData(SS).p_staticCVmin(meetsConstraints),'k','Marker','x');
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
p(ip,iax) = scatter(Vc_l_meetsCon,1e-3*studyData(SS).P_cElec(meetsConstraints),'k','Marker','x');
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
p(ip,iax) = scatter(Vc_l_meetsCon, ...
                100*(studyData(SS).L_c(meetsConstraints)+studyData(SS).L_LPPL(meetsConstraints)),'k','Marker','x');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Mean Power Loss of Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend('charge pump & pipeline combined');
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% charge pump shaft speed
iax = 5;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = scatter(Vc_l_meetsCon,w_c_meetsCon/2/pi*60,'k','Marker','x');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('speed (rpm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Charge Pump Shaft Speed'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% distibution of volume
iax = 6;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = scatter(Vc_l_meetsCon,X_meetsCon,'k','Marker','x');
ip = ip+1;

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('portion', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Portion of Volume at WEC-driven Pump Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])
