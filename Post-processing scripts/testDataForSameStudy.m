files = dir;
nfiles = size(files,1);
dataFileName = "data_parPTO_accum";

% initialize
j = 5; % choose known good data file
load(files(j).name,'-regexp','^(?!out)\w')
Vtotal_last = Vtotal;
X_last = X;
kv_last = kv;

k = 1;
for j = 1:nfiles
display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,dataFileName)
        load(files(j).name,'-regexp','^(?!out)\w')
    
        test(k) = ~any(~(Vtotal == Vtotal_last)) & ~any(~(X == X_last)) & ~any(~(kv == kv_last));
        Vtotal_last = Vtotal;
        X_last = X;
        kv_last = kv;
        k = k + 1;

    end
end

find(test==0) % report files found to be dissimular to previous file