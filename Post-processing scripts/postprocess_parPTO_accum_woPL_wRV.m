%% Collect data from data files
files = dir;
nfiles = size(files,1);
dataFileName = "data_parPTO_accum";
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,dataFileName)
        load(files(j).name,'-regexp','^(?!out)\w')

        % study variables
        studyVar(SS).Vtotal(iVar) = Vtotal_mesh(iVar);
        studyVar(SS).X(iVar) = X_mesh(iVar);
        studyVar(SS).kv(iVar) = kv_mesh(iVar);

        % metrics
        studyData(SS).q_permMean(iVar) = q_permMean;
        studyData(SS).PP_WEC(iVar) = PP_WEC;
        studyData(SS).PP_wp(iVar) = PP_wp;
        studyData(SS).PP_rv(iVar) = PP_rv;
        studyData(SS).L_rv(iVar) = PP_rv/PP_WEC;
        studyData(SS).PP_hinPRV(iVar) = PP_hinPRV;
        studyData(SS).L_hinPRV(iVar) = PP_hinPRV/PP_WEC;
        studyData(SS).PP_roPRV(iVar) = PP_roPRV;
        studyData(SS).L_roPRV(iVar) = PP_roPRV/PP_WEC;
        studyData(SS).PP_pmLoss(iVar) = PP_pmLoss;
        studyData(SS).L_pmLoss(iVar) = PP_pmLoss/PP_WEC;
        studyData(SS).PP_gen(iVar) = PP_gen;
        studyData(SS).X_gen(iVar) = PP_gen/PP_WEC; % proportion of power converted to electricity
        studyData(SS).dpdt_max(iVar) = dpdt_max;
        studyData(SS).dpdt_97(iVar) = dpdt_97(end);

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
            for j = 1:length(notDone(SS).list)
                iVar = notDone(SS).list(j);
                % study variables
                studyVar(SS).Vtotal(iVar) = Vtotal_mesh(iVar);
                studyVar(SS).X(iVar) = X_mesh(iVar);
                studyVar(SS).kv(iVar) = kv_mesh(iVar);

                studyData(SS).q_permMean(iVar) = NaN;
                studyData(SS).PP_WEC(iVar) = NaN;
                studyData(SS).PP_wp(iVar) = NaN;
                studyData(SS).PP_rv(iVar) = NaN;
                studyData(SS).L_rv(iVar) = NaN;
                studyData(SS).PP_hinPRV(iVar) = NaN;
                studyData(SS).L_hinPRV(iVar) = NaN;
                studyData(SS).PP_roPRV(iVar) = NaN;
                studyData(SS).L_roPRV(iVar) = NaN;
                studyData(SS).PP_pmLoss(iVar) = NaN;
                studyData(SS).L_pmLoss(iVar) = NaN;
                studyData(SS).PP_gen(iVar) = NaN;
                studyData(SS).X_gen(iVar) = NaN;
                studyData(SS).dpdt_max(iVar) = NaN;
                studyData(SS).dpdt_97(iVar) = NaN;

            end
        end

    catch
        % just move on
    end
end

%% Transform data to 3D variable mesh
I = length(kv);
J = length(Vtotal);
K = length(X);


test = 1;
for SS = 1:6
    for i = 1:I
        for j = 1:J
            for k = 1:K
                m = J*K*(i-1) + K*(j-1) + k;

                studyData(SS).q_permMean_3D(i,j,k) = studyData(SS).q_permMean(m);
                studyData(SS).PP_WEC_3D(i,j,k) = studyData(SS).PP_WEC(m);
                studyData(SS).PP_wp_3D(i,j,k) = studyData(SS).PP_wp(m);
                studyData(SS).PP_rv_3D(i,j,k) = studyData(SS).PP_rv(m);
                studyData(SS).L_rv_3D(i,j,k) = studyData(SS).L_rv(m);
                studyData(SS).PP_hPRV_3D(i,j,k) = studyData(SS).PP_hinPRV(m);
                studyData(SS).L_hPRV_3D(i,j,k) = studyData(SS).L_hinPRV(m);
                studyData(SS).PP_roPRV_3D(i,j,k) = studyData(SS).PP_roPRV(m);
                studyData(SS).L_roPRV_3D(i,j,k) = studyData(SS).L_roPRV(m);
                studyData(SS).PP_pmLoss_3D(i,j,k) = studyData(SS).PP_pmLoss(m);
                studyData(SS).L_pmLoss_3D(i,j,k) = studyData(SS).L_pmLoss(m);
                studyData(SS).PP_gen_3D(i,j,k) = studyData(SS).PP_gen(m);
                studyData(SS).X_gen_3D(i,j,k) = studyData(SS).X_gen(m);
                studyData(SS).dpdt_max_3D(i,j,k) = studyData(SS).dpdt_max(m);
                studyData(SS).dpdt_97_3D(i,j,k) = studyData(SS).dpdt_97(m);

                Vtotal_3D(i,j,k) = Vtotal_mesh(m);
                X_3D(i,j,k) = X_mesh(m);
                kv_3D(i,j,k) = kv_mesh(m);
                test = (Vtotal_mesh(m) == Vtotal(j)) ...
                    && (X_mesh(m) == X(k)) ...
                    && (kv_mesh(m) == kv(i)) ...
                    && test;
            end
        end
    end
end
if ~test; error('indexing incorrect'); end
clearvars test


%% Find optimal distribution of accumulator volume for each total 
SS = 2;

V_metric = studyVar(SS).Vtotal;
PP_metric = 100*(studyData(SS).PP_rv ...
                + studyData(SS).PP_hinPRV ...
                + studyData(SS).PP_roPRV) ...
                ./studyData(SS).PP_WEC;
dpdt_ub = par.control.dpdt_ROmax;
maxOr97 = 1;
switch maxOr97
    case 1
        dpdt_metric = studyData(SS).dpdt_max;
        varLegLabel = 'max';
        varLabel = 'Maximum';
    case 2
        dpdt_metric = studyData(SS).dpdt_97;
        varLegLabel = '97';
        varLabel = '97th Percentile';
end

% accumulator volume and pipeline diameter
 % find individuals meeting dpdt constraint
dpdt_ub = par.control.dpdt_ROmax;
meetsConstraints = find(dpdt_metric <= dpdt_ub);

 % find non-dominated individuals from set meeting cavitation constraints
non_dominated = paretoFront2D(V_metric(meetsConstraints),'min', ...
               PP_metric(meetsConstraints),'min');
[~, ii_sort] = sort(V_metric(meetsConstraints(non_dominated)));
iiPareto = meetsConstraints(non_dominated(ii_sort));

V_metric_opt = V_metric(iiPareto);
V_metric_meetsCon = V_metric(meetsConstraints);

X_opt = studyVar(SS).X(iiPareto);
X_meetsCon = studyVar(SS).X(meetsConstraints);

kv_opt = studyVar(SS).kv(iiPareto);
kv_meetsCon = studyVar(SS).kv(meetsConstraints);

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

titleString = ['Performance of High-Pressure Circuit Branch',newline,...
                'as a Function of Installed Accumulator Volume',newline,...
                'Sea State ',num2str(SS)];
sgtitle(titleString,...
'Interpreter','latex','FontSize',fontSize+2,'fontname','Times')

n_plots = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Power loss
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

ip = 1;
p(ip,iax) = semilogy(V_metric_opt,PP_metric(iiPareto),'k','Marker','x');
hold on
ip = ip+1;
p(ip,iax) = semilogy(V_metric_opt,100*studyData(SS).L_rv(iiPareto),'b','Marker','o');
ip = ip+1;
p(ip,iax) = semilogy(V_metric_opt,100*(studyData(SS).L_hinPRV(iiPareto)),'r','Marker','^');
ip = ip+1;
p(ip,iax) = semilogy(V_metric_opt,100*(studyData(SS).L_roPRV(iiPareto)),'g','Marker','square');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Mean Power Loss of Normalized to Mean Power Capture'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend('combined','resistive element','PRV at WEC-driven pump','PRV at RO inlet');
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dpdt
iax = 2;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(V_metric_opt,1e-5*dpdt_metric(iiPareto),'k','Marker','x');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('rate of change (bar/s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = [varLabel,' Rate of Change in Pressure at RO Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

% set(leg, 'Location', 'southoutside')
xLim = xlim;
xlim([0 xLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% distibution of volume
iax = 3;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(V_metric_opt,X_opt,'k','Marker','x');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('portion', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Portion of Volume at RO Inlet'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% valve coefficient
iax = 4;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

hold on

ip = 1;
p(ip,iax) = plot(V_metric_opt,kv_opt*sqrt(1000)*1000,'k','Marker','x');
ip = ip+1;

grid on

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('valve coefficient (L/s/kPa^{1/2})', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')

titleString = ['Valve Coefficient'];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

% xLim = xlim;
% xlim([0 xLim(2)])
% yLim = ylim;
% ylim([0 yLim(2)])
linkaxes(ax,'x')
return

%% Plot average power loss from the ripple control valve as a function of total accumulator volume for distribution (color) and valve coefficient (line type)
 % select indices to plot
plotCase = 1;
switch plotCase
 case 1
     iiK = 1:1:K; % distribution
     iiI = 1:2:I;   % valve coeff.
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
    YaxisVar = 100*studyData(SS).PP_rv_3D./studyData(SS).PP_WEC_3D;
    varStr = 'Ripple Control Valve Losses';
  case 2
    YaxisVar = 100*(studyData(SS).PP_roPRV_3D + studyData(SS).PP_hPRV_3D)./studyData(SS).PP_WEC_3D;
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
        ['k_v = ',num2str(kv(iiI(i))*sqrt(1000)*1000,2),'(L/s/kPa^{1/2})']);
end

% plot real data

for k = 1:nK
    for i = 1:nI
        p(k,i) = plot(Vtotal,YaxisVar(iiI(i),:,iiK(k)), ...
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


%% Plot rate of change as a function of total accumulator volume for distribution (color) and valve coefficient (line type)
 % select indices to plot
plotCase = 1;
switch plotCase
 case 1
     iiK = 1:2:K; % distribution
     iiI = 1:3:I;   % valve coeff.
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
    YaxisVar = studyData(SS).dpdt_max_3D;
    varTitle = 'Maximum Rate of Change in Pressure';
  case 2
    YaxisVar = studyData(SS).dpdt_97_3D;
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
        ['k_v = ',num2str(kv(iiI(i))*sqrt(1000)*1000,2),'(L/s/kPa^{1/2})']);
end

% plot real data

for k = 1:nK
    for i = 1:nI
        p(k,i) = plot(Vtotal,1e-3*YaxisVar(iiI(i),:,iiK(k)), ...
            'LineStyle', linestyles{i}, ...
            'Color',color(k,:), ...
            'LineWidth',1);
        p(k,i).HandleVisibility='off';
    end
end

% plot target limit
iLeg = iLeg+1;
plot(Vtotal([1 end]),1e-3*par.control.dpdt_ROmax*[1 1],'--r')
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

%% Plot power loss versus total accumulator volume for pareto optimal 
% designs meeting threshhold on rate of pressure change

XaxisVar = Vtotal_mesh';
YaxisVar = 100*PP_rv_array./PP_WEC_array;
bounds = [1.001 1 0.90 0.75 0.5 0];
N = numel(bounds)-1; % Number of bins with bounds
dpdt_ub = par.control.dpdt_ROmax*bounds;
maxOr97 = 1;
switch maxOr97
    case 1
        dpdt_metric = dpdt_max_array;
        varLegLabel = 'max';
        varLabel = 'Maximum';
    case 2
        dpdt_metric = dpdt_97_array;
        varLegLabel = '97';
        varLabel = '97th Percentile';
end



iDomStart = ones(1,N+1);
iDom = [];
% Loop through bounds
for k = 1:N
    % Find individuals meeting dpdt bounds
    ub = dpdt_ub(k);
    lb = dpdt_ub(k+1);
    [~,meetsConstraints] = find(dpdt_metric <= ub);

    % load objectives
    obj1 = -XaxisVar(meetsConstraints);
    obj2 = -YaxisVar(meetsConstraints);

    % Number of individuals
    n = numel(obj1);

    % Initialize dominance matrix
    dominance = false(n);
    
    % figure
    % scatter(obj1,obj2)
    % hold on
    % Find and mark dominating individuals meeting dpdt criterion
    for i = 1:n
        for j = 1:n
            if obj1(i) >= obj1(j) && obj2(i) >= obj2(j) ...
               && (obj1(i) > obj1(j) || obj2(i) > obj2(j))
                dominance(i,j) = true;
            end
        end
    end
    
    % Identify non-dominated individuals
    non_dominated = find(sum(dominance, 1) == 0)

    % scatter(obj1(non_dominated),obj2(non_dominated))

    % Sort non-dominated individuals based on objective values
    [~, sort_idx] = sort(-obj1(non_dominated));
    Ndom = numel(non_dominated);

    iDomStart(k+1) = iDomStart(k) + Ndom;
    iDom = [iDom meetsConstraints(non_dominated(sort_idx))];
end
iDomStart(k+1) = numel(iDom)+1;

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
iax = 1;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;


iLeg = 0;
for i = 1:N
    iVar = iDom(iDomStart(i):iDomStart(i+1)-1);
    p(iax,i) = plot(XaxisVar(iVar),YaxisVar(iVar));
    p(iax,i).Color = [color(i,:)];
    p(iax,i).Marker = markerType(i);
    p(iax,i).MarkerSize = 5;
    p(iax,i).LineWidth = 1.5;
    hold on

    iLeg = iLeg+1;
    legLabels(iLeg) = convertCharsToStrings( ...
        ['dp/dt|_{',varLegLabel,'} <= ',num2str(1e-3*dpdt_ub(i),3),' kPa/s']);
end

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('power loss (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
switch par.ERUconfig.outlet
    case 1
        ERUstr = '';
    case 0
        ERUstr = 'and ERU Feed Outlet Downstream of Valve';
end
switch par.rvConfig.active
    case 0
        rvStr = 'Passive Ripple Control';
    case 1
        rvStr = 'Active Ripple Control';
end

titleString = ['Power Loss Normalized to Mean Power Capture With ',rvStr,ERUstr,':',newline,...
            varLabel,' Rate of Pressure Change Compared to Limit',newline,...
            'Sea State ',num2str(SS)];
title(titleString,...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')

leg = legend(legLabels);
leg.FontSize = fontSize-1;
leg.FontName = 'Times';
rect = [0.5, -0.2, 0.25, 0.15];
% set(leg, 'Position', rect)
set(leg, 'Location', 'best')
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

for i = 1:N
    iVar = iDom(iDomStart(i):iDomStart(i+1)-1);
    p(iax,i) = plot(XaxisVar(iVar),100*X_mesh(iVar));
    p(iax,i).Color = [color(i,:)];
    p(iax,i).Marker = markerType(i);
    p(iax,i).MarkerSize = 5;
    p(iax,i).LineWidth = 1.5;
    hold on

end

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('portion of volume at RO inlet (x100\%)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distribution of accumulators versus total volume
iax = 3;
ax(iax) = subplot(n_plots,1,iax);
ax(iax).FontName = 'times';
ax(iax).FontSize = fontSize-1;

for i = 1:N
    iVar = iDom(iDomStart(i):iDomStart(i+1)-1);
    p(iax,i) = plot(XaxisVar(iVar),sqrt(1e3)*1000*kv_mesh(iVar));
    p(iax,i).Color = [color(i,:)];
    p(iax,i).Marker = markerType(i);
    p(iax,i).MarkerSize = 5;
    p(iax,i).LineWidth = 1.5;
    hold on

end

xlabel('volume (1000L) ', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylabel('valve coefficient (L/s/kPa^{1/2})', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
xLim = xlim;
xlim([0 xLim(2)])
yLim = ylim;
ylim([0 yLim(2)])

