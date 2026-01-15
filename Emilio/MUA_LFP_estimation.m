exp_paths = ["Z:\PainData\Corrected_Channel_Map\L6\Cortex\20.8.21\KS2", ...
"Z:\PainData\Corrected_Channel_Map\L6\Cortex\26.8.21", ...
"Z:\PainData\Corrected_Channel_Map\Dual\L6_VPL\S1\Nblocks0\th10_2\AUC0pt7\Lambda10", ...
"Z:\PainData\Corrected_Channel_Map\L6\Cortex\m8", ...
"Z:\PainData\Corrected_Channel_Map\L6\Cortex\27.8.21"];
anaesthesia_states = cell(size(exp_paths));
%%
for ce = 1:numel(exp_paths)
    data_dir = exp_paths(ce);
    spike_file = dir(fullfile(data_dir,"*All_all_channels_clean.mat"));
    if numel(spike_file)~=1 
        spike_file = dir(fullfile(data_dir, "*all_channels_clean.mat"));
    end
    load(fullfile(spike_file.folder, spike_file.name));

    unit_flag = sum([sortedData{:,3}] == [1;2]) > 0;
    % unit_flag = [sortedData{:,3}] == 1;
    spike_times = sortedData(unit_flag,2);
    all_spikes = cat(1, spike_times{:});

    bin_files = dir(fullfile(data_dir,'*.bin'));
    if numel(bin_files)>=1
        bin_files = bin_files([bin_files.bytes]==max([bin_files.bytes]));
        exp_duration = (bin_files.bytes/(64*2))/fs;
    else
        exp_duration = max(all_spikes)*1.001;
    end
    fr_per_unit = cellfun(@(x) numel(x)/exp_duration, spike_times);
    %%
    bin_size = 0.1;
    bins = (bin_size/2):bin_size:exp_duration;
    centres = mean([bins(1:end-1);bins(2:end)]);

    muah = histcounts(all_spikes, bins);
    [b1,a1] = butter(3,[0.8,1.8]*2*bin_size,'bandpass');
    [blow, alow] = butter(2,[0.0005, 0.01]*2*bin_size,'bandpass');
    bouts = filtfilt(b1, a1, muah);
    pad_size = ceil(1000/bin_size);
    extra_bout = padarray(bouts(:), pad_size,"symmetric","pre");
    brain_state = filtfilt(blow, alow, abs(extra_bout));
    brain_state(1:pad_size) = [];
    anaesthesia_states{ce} = brain_state;
    %%
    f = figure('PaperSize', [21, 14.8], 'Units', 'centimeters', ...
        'Position', [2 2 21 14.8]);
    t = createtiles(f,1,5);
    axs(1) = nexttile(t,1,[1,4]);
    plot(axs(1),centres, muah, 'Color', (2/3)*ones(1,3), ...
        'DisplayName', 'MUA histogram');
    ylabel(axs(1),'MUA counts')
    yyaxis(axs(1), "right");
    plot(axs(1), centres, zscore(brain_state), 'k','LineWidth', 2, ...
        'DisplayName', 'Anaesthesia state');
    ylabel(axs(1), 'Anaesthesia state')
    legend(axs(1), 'Box', 'off', 'Color', 'none', 'Location', 'best', ...
        'AutoUpdate', 'off')
    yline(axs(1), 0.85, 'k--', 'LineWidth', 1)
    xlabel(axs(1), 'Time [s]'); axs(1).YAxis(2).Color=0.15*ones(1,3);
    xlim(axs(1),[0,exp_duration])

    axs(2) = nexttile(t);
    jit_width = 0.33;
    sc = scatter(axs(2), rand(numel(spike_times),1) * jit_width + 1-(jit_width/2), ...
        fr_per_unit, 'k.', 'MarkerEdgeAlpha',0.75, 'displayname', 'Unit fr');
    hold(axs(2),"on")
    boxchart(axs(2),ones(size(fr_per_unit)),fr_per_unit,"Notch","on", ...
        "BoxEdgeColor","k","BoxFaceColor","none","MarkerStyle","none")
    ln = plot(axs(2),[1,1]+[-1, 1]*jit_width, mean(fr_per_unit)*[1,1], "b", ...
        "LineWidth", 1.5, "DisplayName", ...
        sprintf("Mean %.2g Hz", mean(fr_per_unit)));
    ylabel(axs(2), 'Firing rate [Hz]')
    legend(axs(2), [sc,ln],'Box', 'off', 'Color', 'none', ...
        'Location', 'best', 'AutoUpdate', 'off')
    set(get(axs(2),"XAxis"),"Visible", "off")
    set(axs,"TickDir","out")
    cleanAxis(axs);
    title(t,data_dir,"interpreter","none")
    %%
    saveFigure(f, fullfile(data_dir,'Anaesthesia state estimation and fr'), true, true)
    close(f)
    clearvars -except exp_paths ce anaesthesia_states
end

%%
nas = cellfun(@(x) zscore(x), anaesthesia_states, 'UniformOutput', false);
ths = -3.5:0.01:3.5;
props = zeros(numel(ths), numel(nas));
ii = 1;
for cth = ths
    prop = cellfun(@(x) sum(x>cth)/numel(x), nas);
    props(ii,:) = prop;
    ii = ii + 1;
end