%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_checkValves")
        load(files(j).name,'-regexp','^(?!out)\w')

        eff_wecPump_array(iVar,SS) = eff_wecPump;
        p_min_wp_array(iVar,SS) = p_min_wp;

    end

end

%% Find indices for missing data files

files = dir;
nfiles = size(files,1);
Done = [];
notDone = 1:50*5;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_checkValves")
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
    for j = 1:length(notDone)
        iVar = notDone(j);
        eff_wecPump_array(iVar) = nan;
        p_min_wp_array(iVar) = nan;

    end
    end

catch
    % just move on
end

%% Plot Maximum Rate of Change in Pressure as a function of total accumulator volume
SS = 2;

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

xlabel('flow coefficient, low-pressure (L/s/kPa^{1/2})', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Check Valve Sizing Study Results: Sea State ',num2str(SS)],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')
yyaxis left
semilogx(X*kv*1000*sqrt(1000),eff_wecPump_array(:,SS))
hold on
ylabel('WEC-driven pump efficiency', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
ylim([0 1])

yyaxis right
hold on
semilogx(X*kv*1000*sqrt(1000),1e-3*p_min_wp_array(:,SS))
ylabel('minimum pressure in pump (kPA)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
