%%
close all
figure
idx = ~isnan(HispeedTrials.Base_Ang(92,1).right.whisker1);
time = HispeedTrials.Timestamps{92,1} - HispeedTrials.Timestamps{92,1}(1);
time = time./1000;
plot(time(idx),HispeedTrials.Base_Ang(92,1).right.whisker1(idx))
xlim([0 inf])
ylim([30 120])
axis tight
xlabel('Time since Beam Break [s]')
ylabel('Whisker Base Angle [deg]')
% xlim([0 numel(HispeedTrials.Timestamps{61,1})])
%%
close all
figure
idx = ~isnan(HispeedTrials.Whisk_Curv(64,1).right.whisker2(:,4));
time = HispeedTrials.Timestamps{64,1} - HispeedTrials.Timestamps{64,1}(1);
time = time./1000;
plot(time(idx),HispeedTrials.Whisk_Curv(64,1).right.whisker2(idx,4))
xlim([0 inf])
ylim([30 120])
axis tight
xlabel('Time since Beam Break [s]')
ylabel('Whisker Curvature [1/r]')
% xlim([0 numel(HispeedTrials.Timestamps{61,1})])

%%
close all
count = 1;
for i = 1:size(HispeedTrials,1)
    if sum(~isnan(HispeedTrials.Base_Ang(i,1).right.whisker1)) > 80
        figure, plot(HispeedTrials.Base_Ang(i,1).right.whisker1)
        fprintf("%d: %s\n",count,HispeedTrials.VideoPath(i))
        count = count+1;
    end
end
%%
close all
count = 1;
for i = 1:size(HispeedTrials,1)
    if sum(~isnan(HispeedTrials.Whisk_Curv(i,1).right.whisker1(:,4))) > 100
        figure, plot(HispeedTrials.Whisk_Curv(i,1).right.whisker1(:,4))
        fprintf("%d: %s\n",count,HispeedTrials.VideoPath(i))
        count = count+1;
    end
end

%% Check Excel Sheet for unique cells

excel_sheet = "Z:\Filippo\Animals\Cohort10_23-26\#24\2021-05-03_P3.5_Ruleswitch_EventList.xlsx";

[~,events,RAW]=xlsread(excel_sheet);
RAW(:,4) = cellfun(@(x) strtrim(x),RAW(:,4),'UniformOutput',false);
aa = unique(RAW(:,4));

new_file = "Z:\Filippo\Animals\Cohort10_23-26\#24\2021-05-03_P3.5_Ruleswitch_EventList_Corrected.xlsx";
xlswrite(new_file,RAW)

timestamps = cell2mat(RAW(:,2));
smaller = [];
for i = 2:numel(timestamps)
    if timestamps(i) < timestamps(i-1)
        smaller = [smaller, i];
    end
end

%% Plot mean waveforms
% You need the clWaveforms from the DE_Salience script for that

% figure
% hold on
% for i = 1:size(clWaveforms,1)
for i = 1:20
    figure
    mean_wave = mean(clWaveforms{i,2},2);
    plot(zscore(mean_wave))
end
% hold off

%% Extract mean waveforms and plot the spike width distribution
% You need the clWaveforms from the DE_Salience script for that

widths = nan(1,size(clWaveforms,1));
peak2troughs_amp = nan(1,size(clWaveforms,1));
troughs_loc = nan(1,size(clWaveforms,1));
peak2troughs_dur = nan(1,size(clWaveforms,1));
for unit = 1:size(clWaveforms,1)
    mean_wave = mean(clWaveforms{unit,2},2);
    [~,troughs_loc(unit),widths(unit),peak2troughs_amp(unit)] = findpeaks(-mean_wave,'SortStr','descend','NPeaks',1);
    [~,peak_loc,~,~] = findpeaks(mean_wave);
    if ~isempty(find(peak_loc > troughs_loc(unit),1))
        peak_loc = peak_loc(find(peak_loc > troughs_loc(unit),1));
        peak2troughs_dur(unit) = abs(diff([peak_loc,troughs_loc(unit)]));
    else
        try
            peak_loc = peak_loc(find(peak_loc < troughs_loc(unit),1,'last'));
            peak2troughs_dur(unit) = abs(diff([peak_loc,troughs_loc(unit)]));
        catch
            peak2troughs_dur(unit) = nan;
        end
    end
end

% Widths are given in data points. Calculate in usecs
framerate = fs; % in Hz
widths_usec = (widths/framerate)*1000000;
peak2troughs_usec = (peak2troughs_dur/framerate)*1000000;
figure
hold on
% h(1) is the handle to the histogram, and h(2) is the handle to the density curve
h = histfit(widths_usec,50,'kernel');
localmin = h(2).XData(islocalmin(h(2).YData));
localmin = localmin(1);
ylims = ylim;
plot([localmin, localmin], [ylims(1),ylims(2)], 'r:', 'LineWidth', 3)
text(localmin, 0.9*ylims(2),sprintf('%.2f µm  ',localmin),'HorizontalAlignment','right')
set(gca, 'XLimSpec', 'Tight')
xlabel('Full Width at Half Maximum [µs]')
ylabel('Cell count')
title('Distribution of waveform widths across all units')
hold off

% Gives the total amount of units with a mean waveform smaller than the
% local min (i.e., putative interneurons)
put_in = sum(h(1).YEndPoints(1:find(h(1).XData < localmin,1,'last')));
% Gives the total amount of units with a mean waveform greater than the
% local min (i.e., excitatory interneurons)
put_ex = sum(h(1).YEndPoints(find(h(1).XData < localmin,1,'last')+1:end));

figure
pie([put_in,put_ex],'%.1f%%')

figure
s_1 = scatterhist(widths_usec,peak2troughs_usec,'NBins',30,'Kernel','overlay');
axis padded

% Calculate principal component
coeff = pca([widths_usec',peak2troughs_usec']);

ylims = ylim;
xvals = linspace(coeff(1)+coeff(2)*ylims(1),coeff(1)+coeff(2)*ylims(2));
hold on, plot(xvals,coeff(1)+coeff(2)*xvals)

xlabel('Full Width at Half Maximum [µs]')
ylabel('Trough to Second Peak [µs]')
title('Scatter of Waveform Characterstics')

% Rotate along the first principal component
theta = atan(coeff(2));
theta_deg = (theta/pi)*180;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

v = [widths_usec' peak2troughs_usec'];
vr = v*R;
x = vr(:,1);
y = vr(:,2);

% figure
% s_2 = scatterhist(x,y,'NBins',30);

% % Plot PC line as reference
% v_pc = [linspace(0,800)' (coeff(1)+coeff(2)*linspace(0,800))'];
% vr_pc = v_pc*R;
% x_pc = vr_pc(:,1);
% y_pc = vr_pc(:,2);
% hold on, plot(x_pc,y_pc), hold off

%%
figure

subplot(2,1,1)
h2 = histfit(x,30,'kernel');

localmin = h2(2).XData(islocalmin(h2(2).YData));
ylims = ylim;
try 
    hold on, plot([localmin, localmin], [ylims(1),ylims(2)], 'r:', 'LineWidth', 3)
catch
end

xticklabels({})
ylabel('Cell Count')
axis tight
title('Cell Distribution along the PC1')

subplot(2,1,2)
s_3 = scatter(widths_usec,peak2troughs_usec);
axis padded
ylims = ylim;
xvals = linspace(coeff(1)+coeff(2)*ylims(1),coeff(1)+coeff(2)*ylims(2));
hold on, plot(xvals,coeff(1)+coeff(2)*xvals)

title('Scatter of Waveform Characterstics with the PC1')
xlabel('Full Width at Half Maximum [µs]')
ylabel('Trough to Second Peak [µs]')

%%
figure
s_3 = scatter(widths_usec,peak2troughs_usec);

a1 = gca;

a3 = axes('position',get(a1,'position'));
set(a3,'color','none','ycolor','none','xcolor','none');
% a3.Visible = 'off';
title('Scatter of Waveform Characterstics rotated along PC1')

a2 = axes('position',get(a1,'position'));
h2 = histogram(flip(x),'NumBins',30);
pos1 = get(a1,'Position');
posfig = get(gcf,'Position');

rlim=[min([x(:);y(:)]),max([x(:);y(:)])];rlim(2)=rlim(2)+rlim(2)-rlim(1);

% tmp = min(pos1([3,4]).*posfig([3,4]));
% pos1(3) = tmp/posfig(3);
% pos1(4) = tmp/posfig(4);
% pos1(3:4) = [0.7 0.7];
pos1 = [0.23 0.11 0.6 0.4];
pos2 = [0.23 0.11 0.6 0.2];
% pos2 = [pos1(1)+(pos1(3))/2-pos1(3)/6, pos1(2)+(pos1(4))/2-pos1(4)/6, 0.6, 0.4];

set(a1,'ActivePositionProperty','position','Position',pos1,'xlim',rlim,...
    'ylim',rlim,'color','none');
set(a2,'ActivePositionProperty','position','Position',pos2,'color','none');
set(a2,'CameraViewAngle',0,'View',[theta_deg-45,90],'box','off','ycolor','none',...
    'xtick',[],'Nextplot','add');
set(gcf,'CurrentAxes',a1)

xlabel('Full Width at Half Maximum [µs]')
ylabel('Peak to Trough [µs]')

% title('Scatter of Waveform Characterstics rotated along PC1')
axis padded

%%
modelfun = @(b,x) b(1) + b(2) * x + b(3) * exp(-(x(:, 1) - b(4)).^2/b(5)) + b(6) * exp(-(x(:, 1) - b(7)).^2/b(8));  
beta0 = [6, 0.1, 35, 10, 3, 25, 50, 9]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);
% Now the model creation is done and the coefficients have been determined.
% YAY!!!!
% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl.Coefficients{:, 'Estimate'};
% Let's do a fit, but let's get more points on the fit, beyond just the widely spaced training points,
% so that we'll get a much smoother curve.
X = linspace(min(X), max(X), 1920); % Let's use 1920 points, which will fit across an HDTV screen about one sample per pixel.
% Create smoothed/regressed data using the model:
yFitted = coefficients(1) + coefficients(2) * X + coefficients(3) * exp(-(X - coefficients(4)).^2 / coefficients(5)) + coefficients(6) * exp(-(X - coefficients(7)).^2 / coefficients(8));



