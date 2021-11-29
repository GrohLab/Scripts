function Classifiers = getWaveformClassifiers(mean_wf, fs)


[Nt, Ncl] = size(mean_wf);
% Window multiplication to ensure we look at the center part of the
% waveform
%  winMult = hann(Nt)';
%  mean_wf = mean_wf .* winMult';

tx = (0:Nt-1)'./fs - Nt/(2*fs);
tcp = getWaveformCriticalPoints(mean_wf, fs);
tcp = cellfun(@(x) x - Nt/(2*fs), tcp, 'UniformOutput', 0);
b50 = min(mean_wf) + range(mean_wf)./2;
ampCp = cell(Ncl,1);
tpd = zeros(Ncl,1); bhad = tpd;
frstPk = zeros(Ncl,1);
secPk = zeros(Ncl,1);

for ccl = 1:Ncl
    ampCp(ccl) = {interp1(tx, mean_wf(:,ccl), tcp{ccl,1}, 'spline')};
    [~, troughSub] = max(abs(ampCp{ccl}));
    [sP, secPeakSub] = max(abs(ampCp{ccl}(troughSub+1:end)));
    if isempty(secPeakSub)
        % We have to think of a math solution for the so called 'stable
        % region' or 'period'.
        fprintf(1, 'Cluster %d is a bitch\n', ccl)
        figure(ccl); plot(tx, mean_wf(:,ccl))
        continue;
    end
    fP = max(ampCp{ccl}(1:troughSub-1));
    tpd(ccl) = diff(tcp{ccl,1}([troughSub, troughSub+secPeakSub]));
    if ampCp{ccl}(troughSub) > 0
        halfFlags = mean_wf(:,ccl) >= b50(ccl);
    else
        halfFlags = mean_wf(:,ccl) <= b50(ccl);
    end
    frstSub = find(halfFlags, 1, 'first');
    if frstSub <= 1
        Classifiers = NaN(1,4);
    else
        
        lstSub = find(halfFlags, 1, 'last');
        mdl = fit_poly(tx([frstSub-1, frstSub]), mean_wf([frstSub-1, frstSub],ccl),...
            1); frstXg = (b50(ccl) - mdl(2))/mdl(1);
        mdl = fit_poly(tx([lstSub-1, lstSub]), mean_wf([lstSub-1, lstSub],ccl),...
            1); lstXg = (b50(ccl) - mdl(2))/mdl(1);
        bhad(ccl) = lstXg - frstXg;
    end
    if isempty(fP)
        frstPk(ccl) = 0;
    else
        frstPk(ccl) = fP;
    end
    
    
    secPk(ccl) = sP;
    
    
end
Classifiers = [tpd, bhad, frstPk, secPk];
end