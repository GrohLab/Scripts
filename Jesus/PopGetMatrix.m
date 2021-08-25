ClustersL1OnsetVPM = nan(2000,size(PoPrelativeSpkTmsL1,1));
for i = 1:size(PoPrelativeSpkTmsL1(:,1),1)
    Spikes=cell2mat(PoPrelativeSpkTmsL1(i,1:end));
    N= size(Spikes,2);
    ClustersL1OnsetVPM (1:N,i) = Spikes;
end

ClustersL10OnsetVPM = nan(2000,size(PoPrelativeSpkTmsL10,1));
for i = 1:size(PoPrelativeSpkTmsL10(:,1),1)
    Spikes=cell2mat(PoPrelativeSpkTmsL10(i,1:end));
    N= size(Spikes,2);
    ClustersL1OnsetVPM (1:N,i) = Spikes;
end

ClustersL50OnsetVPM = nan(2000,size(PoPrelativeSpkTmsL50,1));
for i = 1:size(PoPrelativeSpkTmsL50(:,1),1)
    Spikes=cell2mat(PoPrelativeSpkTmsL50(i,1:end));
    N= size(Spikes,2);
    ClustersL50OnsetVPM (1:N,i) = Spikes;
end

ClustersL100OnsetVPM = nan(2000,size(PoPrelativeSpkTmsL100,1));
for i = 1:size(PoPrelativeSpkTmsL100(:,1),1)
    Spikes=cell2mat(PoPrelativeSpkTmsL100(i,1:end));
    N= size(Spikes,2);
    ClustersL100OnsetVPM (1:N,i) = Spikes;
end

ClustersL200OnsetVPM = nan(2000,size(PoPrelativeSpkTmsL200,1));
for i = 1:size(PoPrelativeSpkTmsL200(:,1),1)
    Spikes=cell2mat(PoPrelativeSpkTmsL200(i,1:end));
    N= size(Spikes,2);
    ClustersL200OnsetVPM (1:N,i) = Spikes;
end

ClustersControlOnsetVPM = nan(2000,size(PoPrelativeSpkTmsControl,1));
for i = 1:size(PoPrelativeSpkTmsControl(:,1),1)
    Spikes=cell2mat(PoPrelativeSpkTmsControl(i,1:end));
    N= size(Spikes,2);
    ClustersControlOnsetVPM (1:N,i) = Spikes;
end