fnOpts = {'UniformOutput', false};
Ns_puc = arrayfun(@(c) ...
    arrayfun(@(u) ...
    cellfun(@(t) all([~isempty(t), t>0, t<0.05]), ...
    firstSpkStruct(c).FirstSpikeTimes(u,:)), ...
    1:size(firstSpkStruct(c).FirstSpikeTimes, 1), fnOpts{:}), ...
    1:length(firstSpkStruct), fnOpts{:});

fs_puc = arrayfun(@(c) ...
    arrayfun(@(u) ...
    [repmat(c, sum(Ns_puc{c}{u}), 1), ...
    repmat(u, sum(Ns_puc{c}{u}), 1), ...
    [firstSpkStruct(c).FirstSpikeTimes{u,Ns_puc{c}{u}}]'], ...
    1:size(firstSpkStruct(c).FirstSpikeTimes, 1), fnOpts{:}), ...
    1:length(firstSpkStruct), fnOpts{:});

fs_puc = cellfun(@(c) cat(1, c{ cellfun(@(c2) ~isempty(c2),c) }), fs_puc, fnOpts{:});
fs_puc = cat(1, fs_puc{:});