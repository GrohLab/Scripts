function plotTriggeredDoseResponse(ID, sortedData, condIndices, Conditions, fs)



            % atm only works on block conditions (plotUnitInfo omits these for now)....come back to this

% Naming the figure
if class(ID) == 'cell'
    figureName = ID{1};
    lngth = 1; % generalise this!
else
    figureName = ID;
    lngth = 1;
end

% Finding the different frequencies - probably lots of redundancy
rates = zeros(lngth, length(condIndices));
sds = zeros(lngth, length(condIndices));
timeBin = 5;
conds = zeros(1,length(condIndices));

alltext = cat(2, Conditions.name);
undscore = strfind(alltext, "_");
allHzind = strfind(alltext, "Hz");
freqConds = zeros(length(allHzind),1);
hzUnd = zeros(length(allHzind),1);
for a = 1:length(allHzind)
    hzUnd(a) = max(undscore(undscore<allHzind(a)));
    freqConds(a) = str2double(alltext(hzUnd(a)+1:allHzind(a)-1));
end
freqConds = unique(freqConds);
for a = 1:length(freqConds)
    FreqConds{a} = [num2str(freqConds(a)), 'Hz'];
end

% Getting the mean rates and SDs for the considered Conditions 
for b = 1:lngth
    spkInd = ismember(sortedData(:,1), ID);
    train = sortedData(spkInd, 2);
    Cond = cell(1, length(condIndices));
    for a = 1:length(condIndices)
        Cond{a} = Conditions(condIndices(a)).Triggers;
        rate = TriggeredRates(train, Cond{a},fs, timeBin);
        rates(b,a) = mean(rate{:});
        sds(b,a) = std(rate{:});
        for i = 1:length(FreqConds)
            
            if contains(Conditions(condIndices(a)).name, FreqConds{i}) % sorting the rates into frequency groups
                conds(a) = i;
            end
        end
        if contains(Conditions(condIndices(a)).name, ["pulse", "control"], 'IgnoreCase', true)
            conds(a) = i+1; % looking for continuous pulse conditions
        end
    end
end
condInd = unique(conds);
condFlg = false(length(condInd), length(conds));
for a = 1:length(condInd)
    condFlg(a,:) = conds == condInd(a);
end
consConds = sum(condFlg,2) ~= false;


cc = find(consConds);
legText = [];
legendCell = cell(length(condInd),1);
for a = 1:length(condInd)-1
    legendCell{a} = [FreqConds{a}, legText];
end
legendCell{a+1} = ['Pulse', legText];
ax = gca;


% Sorting by Laser Intensity

xnames = cell(1, length(condIndices));
pwrs = zeros(1, length(condIndices));
for a = 1:length(condIndices)
    mWind = strfind(Conditions(condIndices(a)).name, 'mW');
    undInds = strfind(Conditions(condIndices(a)).name, '_');
    mwUnd = undInds(undInds<mWind);
    mWundInd = mwUnd(end);
    pwrs(a) = str2double(Conditions(condIndices(a)).name(mWundInd+1:mWind-1));
    
    Conditions(condIndices(a)).name(undInds) = ' ';
    Conditions(condIndices(a)).name = Conditions(condIndices(a)).name(undInds(1)+1:undInds(end)-1);
    if strfind(Conditions(condIndices(a)).name, 'sec')
        Conditions(condIndices(a)).name = Conditions(condIndices(a)).name(6:end);
        Conditions(condIndices(a)).name(1) = 'P';
    end
    xnames{a} = Conditions(condIndices(a)).name;
    
    
end
pwrs = unique(pwrs);
pwrCell = cell(length(pwrs),1);
for a = 1:length(pwrs)
    pwrCell{a} = num2str(pwrs(a));
end


% Plotting

hold on
x = 1:length(pwrs);
offset = 0.1;
for a = 1:length(cc)
    errorbar(ax, x + offset*(a-1), rates(condFlg(cc(a),:)), sds(condFlg(cc(a),:)), 'vertical', 'LineWidth', 3);
end

% Formatting

legend(legendCell, 'Location', 'northwest')
xlabel('Laser Intensity _{(mW)}');
ax.XLim(1) = 0;
ax.XLim(2) = length(pwrs) + 1;
ax.XTick = [1:1:length(pwrs)];
ax.XTickLabel = pwrCell;
ax.YLim(1) = 0;
ylabel('Firing Rate _{(Hz)}');
ax.FontSize = 20;
ax.FontName = 'Arial';
end