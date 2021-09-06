% Unit Depth Re-ordering besed on corrected channel maps.
depths = [0:20:1275]';
for ch = 1:length(chanMap0ind)
    chNo = chanMap0ind(ch);
    ind = find(rhdInput == chNo);
    clInfo.depth(clInfo.ch == chNo) = depths(ind);
end