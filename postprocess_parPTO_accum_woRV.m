%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_accum_woRV")
        load(files(j).name,'-regexp','^(?!out)\w')
        for iX = 1:nVar3
            % study variables
            studyVar(SS).Vc_h(iVar,iX) = Vc_h_mesh(iVar);
            studyVar(SS).X(iVar,iX) = X(iX);
            studyVar(SS).d_line(iVar,iX) = d_HPPL_mesh(iVar);

            % metrics
            studyData(SS).q_permMean(iVar,iX) = q_permMean(iX);
            studyData(SS).PP_WEC(iVar,iX) = PP_WEC(iX);
            studyData(SS).PP_wp(iVar,iX) = PP_wp(iX);
            studyData(SS).PP_rv(iVar,iX) = PP_rv(iX);
            studyData(SS).PP_hinPRV(iVar,iX) = PP_hinPRV(iX);
            studyData(SS).PP_roPRV(iVar,iX) = PP_roPRV(iX);
            studyData(SS).dpdt_max(iVar,iX) = dpdt_max(iX);
            studyData(SS).dpdt_97(iVar,iX) = dpdt_97(iX);
            studyData(SS).P_HPPL(iVar,iX) = P_HPPL(iX);
            studyData(SS).L_HPPL(iVar,iX) = L_HPPL(iX);
        end
    end
end

%% Find indices for missing data files

files = dir;
nfiles = size(files,1);
Done = [];
notDone = 1:nVar;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_accum")
        load(files(j).name,'-regexp','^(?!out)\w')
        [r,c,val] = find(notDone==iVar);
        notDone = [notDone(1:c-1), notDone(c+1:end)];
        Done = [Done, iVar];
    end
end

try 
    doneArrayStr = num2str(Done(1));
    for j = 2:length(Done)
        doneArrayStr = append(arrayStr,[',',num2str(Done(j))]);
    end
catch
    % just move on
end

try
    jobArrayStr = num2str(notDone(1));
    for j = 2:length(notDone)
        jobArrayStr = append(jobArrayStr,[',',num2str(notDone(j))]);
    end
    
    if 1
        nanArray = nan*ones(numel(X),1);
        for j = 1:length(notDone)
            iVar = notDone(j);
            studyData.q_permMean(iVar,:) = nanArray;
            studyData.PP_WEC(iVar,:) = nanArray;
            studyData.PP_wp(iVar,:) = nanArray;
            studyData.PP_rv(iVar,:) = nanArray;
            studyData.PP_hinPRV(iVar,:) = nanArray;
            studyData.PP_roPRV(iVar,:) = nanArray;
            studyData.dpdt_max(iVar,:) = nanArray;
            studyData.dpdt_97(iVar,:) = nanArray;
            studyData.P_HPPL(iVar,:) = nanArray;
            studyData.L_HPPL(iVar,:) = nanArray;
        end
    end

catch
    % just move on
end

%% Find optimal distribution of accumulator volume for each total 
SS = 2;

V_metric = studyVar(SS).Vc_h(:);
PP_metric = 100*(studyData(SS).PP_rv(:) ...
                + studyData(SS).PP_roPRV(:) ...
                + studyData(SS).PP_hinPRV(:) ...
                + studyData(SS).P_HPPL(:)) ...
                ./studyData(SS).PP_WEC(:);
dpdt_ub = 0.7e5; % [Pa/s]
maxOr97 = 1;
switch maxOr97
    case 1
        dpdt_metric = studyData(SS).dpdt_max(:);
        varLegLabel = 'max';
        varLabel = 'Maximum';
    case 2
        dpdt_metric = studyData(SS).dpdt_97(:);
        varLegLabel = '97';
        varLabel = '97th Percentile';
end


% Find individuals meeting dpdt bounds
meetsConstraints = find(dpdt_metric <= dpdt_ub);

 % find non-dominated individuals from set meeting dpdt criterion
non_dominated = paretoFront2D(V_metric(meetsConstraints),'min', ...
                              PP_metric(meetsConstraints),'min');
[~, ii_sort] = sort(V_metric(meetsConstraints(non_dominated)));
ii = meetsConstraints(non_dominated(ii_sort));
V_metric_opt = V_metric(ii);
PP_metric_opt = PP_metric(ii);
X_opt = studyVar(SS).X(ii);
d_HPPL_opt = studyVar(SS).d_line(ii);

clearvars meetsConstraints non_dominated ii ii_sort

%% Plot Pereto optimal results
% Plot power loss versus total accumulator volume for pareto optimal
% designs meeting threshhold on rate of pressure change



black = [0 0 0];
maroon = [122 0 25]/256;
gold = [255 204 51]/256;
blue = [0 75 135]/256;
orange = [226 100 55]/256;
green = [63 150 87]/256;
color = [maroon; gold; blue; orange; green];
markerType = '.ox*^s';

bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 9;
lineWidth = 1;

clearvars leg

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

n_plots = 3;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% metrics
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

p(iax) = plot(V_metric_opt,PP_metric_opt);



xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

switch par.rvConfig.active
    case 0
        rvStr = 'Passive Ripple Control';
    case 1
        rvStr = 'Active Ripple Control';
end

titleString = ['Power Loss Normalized to Mean Power Capture With ',rvStr,':',newline,...
            varLabel,' Rate of Pressure Change Compared to Limit',newline,...
            'Sea State ',num2str(SS)];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distribution of accumulators versus total volume
iax = 2;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

p(iax) = plot(V_metric_opt,100*X_opt);

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('portion of volume at RO inlet (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipeline diam. versus total volume
iax = 3;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

p(iax) = plot(V_metric_opt,100*d_HPPL_opt);

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('pipe diameter (cm)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])


%% Transform data to 3D variable mesh
I = length(kv);
J = length(Vc_h);
K = length(X);


test = 1;
for i = 1:I
    for j = 1:J
        for k = 1:K
            m = J*K*(i-1) + K*(j-1) + k;

            q_permMean_3D(i,j,k) = studyData.q_permMean;
            PP_WEC_3D(i,j,k) = studyData.PP_WEC_array(m);
            PP_wp_3D(i,j,k) = studyData.PP_wp_array(m);
            PP_rv_3D(i,j,k) = studyData.PP_rv_array(m);
            PP_hPRV_3D(i,j,k) = studyData.PP_hinPRV_array(m);
            PP_roPRV_3D(i,j,k) = studyData.PP_roPRV_array(m);
            dpdt_max_3D(i,j,k) = studyData.dpdt_max_array(m);
            dpdt_97_3D(i,j,k) = studyData.dpdt_97_array(m);
            P_HPPL_3D(i,j,k) = studyData.P_HPPL(m);
            L_HPPL_3D(i,j,k) = studyData.L_HPPL(m);

            Vc_h_3D(i,j,k) = Vc_h_mesh(m);
            X_3D(i,j,k) = X_mesh(m);
            d_HPPL_3D(i,j,k) = d_HPPL_mesh(m);
            test = (Vc_h_mesh(m) == Vc_h(j)) ...
                && (X_mesh(m) == X(k)) ...
                && (d_HPPL_mesh(m) == d_HPPL(i)) ...
                && test;
        end
    end
end

if ~test; error('indexing incorrect'); end
clearvars test
%% Plot average power loss as a function of total accumulator volume for
% distribution (color) and pipeline ID (line type)
 % select indices to plot
plotCase = 1;
switch plotCase
 case 1
     iiK = 3:1:K-2; % distribution
     iiI = 2:1:I;   % pipeline diam.
 case 2
     iiK = 1:1:K;
     iiI = 4;
 case 3
     iiK = 3;
     iiI = 2:I;
end
nK = length(iiK);
nI = length(iiI);

 % select variable to plot
switch 1
  case 1
    PP_metric = 100*L_HPPL_3D;
    varStr = 'Pipeline Losses';
  case 2
    PP_metric = 100*(PP_roPRV_3D + PP_hinPRV_3D)./PP_WEC_3D;
    varStr = 'PRV Losses';
end

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

% dummy plots for legend
iLeg = 0;

 % color - total volume
for k = 1:nK
    scatter(-99*[1, 0.5],-99*[1, 0.5],50, ...
        'filled','s','LineWidth',2,...
        'MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:));
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['fraction at RO inlet = ',num2str(X(iiK(k)))]);
end

for i = 1:nI
    plot(-99*[1, 0.5],-99*[1, 0.5],'k','LineStyle', linestyles{i});
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['ID = ',num2str(d_HPPL(iiI(i))*1e2,2),'(cm)']);
end

% plot real data

for k = 1:nK
    for i = 1:nI
        p(k,i) = plot(Vc_h,PP_metric(iiI(i),:,iiK(k)), ...
            'LineStyle', linestyles{i}, ...
            'Color',color(k,:), ...
            'LineWidth',1);
        p(k,i).HandleVisibility='off';
    end
end

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

switch par.rvConfig.active
    case 0
        rvStr = 'Passive Ripple Control';
    case 1
        rvStr = 'Active Ripple Control';
end

titleString = ['Mean Power Loss Normalized to Mean Power Capture',newline...
                'With ',rvStr,': ',varStr,newline,...
                'Sea State ',num2str(SS)];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels);
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
% set(leg, 'Location', 'southoutside')
xLim = xlim;
xlim([0 xLim(2)])
% ylim([0 1e-3*max(Y(:))])


%% Plot rate of change as a function of total accumulator volume for
% distribution (color) and pipeline ID (line type)
 % select indices to plot
plotCase = 1;
switch plotCase
 case 1
     iiK = 3:1:K-2; % distribution
     iiI = 2:1:I;   % pipeline diam.
 case 2
     iiK = 1:1:K;
     iiI = 4;
 case 3
     iiK = 3;
     iiI = 2:I;
end
nK = length(iiK);
nI = length(iiI);

 % select variable to plot
maxOr97 = 1;
switch maxOr97
  case 1
    dpdt_metric = dpdt_max_3D;
    varTitle = 'Maximum Rate of Change in Pressure';
  case 2
    dpdt_metric = dpdt_97_3D;
    varTitle = '97th Percentile Rate of Change in Pressure';
end

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

% dummy plots for legend
iLeg = 0;

 % color - total volume
for k = 1:nK
    scatter(-99*[1, 0.5],-99*[1, 0.5],50, ...
        'filled','s','LineWidth',2,...
        'MarkerEdgeColor',color(k,:),'MarkerFaceColor',color(k,:));
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['fraction at RO inlet = ',num2str(X(iiK(k)))]);
end

for i = 1:nI
    plot(-99*[1, 0.5],-99*[1, 0.5],'k','LineStyle', linestyles{i});
    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['ID = ',num2str(d_HPPL(iiI(i))*1e2,2),'(cm)']);
end

% plot real data

for k = 1:nK
    for i = 1:nI
        p(k,i) = plot(Vc_h,1e-3*dpdt_metric(iiI(i),:,iiK(k)), ...
            'LineStyle', linestyles{i}, ...
            'Color',color(k,:), ...
            'LineWidth',1);
        p(k,i).HandleVisibility='off';
    end
end

% plot target limit
iLeg = iLeg+1;
plot(Vc_h([1 end]),1e-3*par.control.dpdt_ROmax*[1 1],'--r')
legLabels(iLeg) = convertCharsToStrings( ...
        ['target limit']);

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('rate of change in pressure (kPa/s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

switch par.rvConfig.active
    case 0
        rvStr = 'Passive Ripple Control';
    case 1
        rvStr = 'Active Ripple Control';
end

titleString = [varTitle,' With ',rvStr,':',newline,...
                'Sea State ',num2str(SS)];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels);
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
% ylim([0 1e-3*1.5*par.control.dpdt_ROmax])