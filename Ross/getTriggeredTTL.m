function[stims, r, csNames] = getTriggeredTTL(cst, delayFlags, trigNames, clr, Nt)
csNames = trigNames;
stims = mean(cst(:,:,delayFlags),3);
    stims = stims - median(stims,2);
    for cs = 1:size(stims,1)
        if abs(log10(var(stims(cs,:),[],2))) < 13
            [m,b] = lineariz(stims(cs,:),1,0);
            stims(cs,:) = m*stims(cs,:) + b;
        else
            stims(cs,:) = zeros(1,Nt);
        end
    end
    [r,c] = size(stims);
    if r < c
        stims = stims';
    end
    
    if sum(clr(3)' & contains(trigNames, "Laser")) ~= false
        csNames = csNames(contains(trigNames, "Laser"));
        r = find(contains(trigNames, "Laser"));
    else
        csNames = csNames(~contains(trigNames, "Laser"));
        r = find(~contains(trigNames, "Laser"));
    end
end