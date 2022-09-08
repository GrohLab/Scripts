%% Trial by trial ephys-roller speed relationship
% PSTH per unit per trial, per condition.
PSTH_unit_trial = getPSTH_perU_perT(relativeSpkTmsStruct, configStructure);

% Median of firing rate in the responsive window (evoked)
respFlag = psthTx(:) > responseWindow;
respFlag = xor(respFlag(:,1), respFlag(:,2));
evoked_median_unit_fr = cellfun(@(x) squeeze(median(x(:, respFlag, :), 2)), ...
    PSTH_unit_trial, fnOpts{:});

% Trial identification between behaviour and ephys
cum_ephys_trial_ID = arrayfun(@(x) cumsum(delayFlags(:,x)), ...
    1:size(delayFlags,2), fnOpts{:});
cum_speed_trial_ID = arrayfun(@(x) cumsum(xdf(:,x)), ...
    1:size(xdf, 2), fnOpts{:});
considered_ephys_trials = arrayfun(@(x) cum_ephys_trial_ID{x}(xdf(:,x)), ...
    1:size(xdf, 2), fnOpts{:});

[mdl, ~, r2] = cellfun(@(x, y, v) ...
    arrayfun(@(u) fit_poly(x(v,u), y, 1), ...
    (1:size(x,2))', fnOpts{:}), ...
    evoked_median_unit_fr(:), mvpt, considered_ephys_trials(:), fnOpts{:});

linear_models = cellfun(@(x) cat(2, x{:}), mdl, fnOpts{:});
r2 = cellfun(@(x) cat(2, x{:}), r2, fnOpts{:});
%{
for ccl = 1:Ncl
    figure('Name', sprintf('%s', gclID{ccl}));
    scatter(evoked_median_unit_fr{1}(considered_ephys_trials{1},ccl), ...
        mvpt{1}, 'k.')
end
%}
%% Population encoding
PSTH_trial = cellfun(@(x) mean(x, 3), PSTH_unit_trial, fnOpts{:});

% population_model
evoked_median_pop_fr = cellfun(@(x, v) median(x(v,respFlag), 2), ...
    PSTH_trial, considered_ephys_trials,  fnOpts{:});

time_color = arrayfun(@(x) jet(sum(xdf(:,x))), 1:size(xdf,2), fnOpts{:});
figure('Name','Population mean','Colormap',time_color{1})
scatter(evoked_median_pop_fr{1}, mvpt{1}, [], time_color{1}, 'filled', ...
    'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.7)
xlabel('Population firing rate [Hz]')
ylabel('Maximum roller speed [cm/s]')
title('Population firing rate vs roller speed')
cbObj = colorbar('Box', 'off'); cbObj.Label.String = 'Normalized time';
set(gca, 'Color', 'none', 'Box', 'off', 'Clipping','off')

%% 
PSTH_unit_trial = getPSTH_perU_perT(relativeSpkTmsStruct, configStructure);
PSTH_trial = cellfun(@(x) mean(x, 3), PSTH_unit_trial, fnOpts{:});