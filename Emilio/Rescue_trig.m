for ct = lSub(lSub(:,1) > 1e3*fs,:)'
    trig(2,ct(1):ct(2)) = trig(2, ct(1):ct(2)) + 2^15 - 230;
end
%%
for ct = (Conditions(2).Triggers(Conditions(2).Triggers(:,1)>expSamples(1),:) - expSamples(1))'
    trig(2,ct(1):ct(2)) = trig(2, ct(1):ct(2)) + int16(2^15 - 230);
end