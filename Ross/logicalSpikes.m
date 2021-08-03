logicalSpks = false(Ncl, Ns);
for cl = 1:Ncl
for spk = 1: length(spkSubs{cl})
    ind = spkSubs{cl}(spk);
logicalSpks(cl,ind) = true;
end
end