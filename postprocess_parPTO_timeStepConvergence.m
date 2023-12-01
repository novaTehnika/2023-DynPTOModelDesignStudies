%% Collect data from data files
files = dir;
nfiles = size(files,1);
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_timeStepConvergence")
        load(files(j).name,'-regexp','^(?!out)\w')

        Ebal_error_array(iVar,SS) = Ebal_error;
        Vbal_error_array(iVar,SS) = Vbal_error;

    end

end

%% Find indices for missing data files

files = dir;
nfiles = size(files,1);
Done = [];
notDone = 1:50*5;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,"data_parPTO_timeStepConvergence")
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

%% Plot normalized energy and volume balance error as a function of timestep
SS = 3;

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

xlabel('time step (s)', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
title(['Normalized Energy and Volume Balance Study: Sea State',num2str(SS)],...
'Interpreter','latex','FontSize',fontSize,'fontname','Times')
yyaxis left
semilogx(MaxStep,abs(Ebal_error_array(:,SS)))
hold on
ylabel('energy balance error, normalized', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')


yyaxis right
hold on
semilogx(MaxStep,abs(Vbal_error_array(:,SS)))
ylabel('volume balance error, normalized', ...
'Interpreter','latex','FontSize',fontSize-1,'fontname','Times')
