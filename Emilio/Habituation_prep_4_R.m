miceArray = [];
for cm = 1:numel(resTable)
    [Np, Ns] = size(resTable{cm}.BehaviourIndices);
    sess_id = repmat(1:Ns, Np, 1); sess_id = sess_id(:);
    mouse_id = cm + zeros(Np*Ns, 1);
    puff_strength = ...
        cellfun(@(f) textscan(f, "%fbars"), ...
        resTable{cm}.Properties.RowNames);
    if any(cellfun(@isempty, puff_strength))
        puff_strength{cellfun(@isempty, puff_strength)} = -1;
    end
    puff_strength = repmat(cat(1, puff_strength{:}), Ns, 1);
    bi_values = resTable{cm}.BehaviourIndices(:);
    miceArray = [miceArray;
        mouse_id, sess_id, puff_strength, bi_values];
end