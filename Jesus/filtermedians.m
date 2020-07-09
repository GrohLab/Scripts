
a=yL200On-yControlOn
BiggerMedianClusters=a>0.000;
SmallerMedianClusters=a<0.000;
SmallerMedianClusters=SmallerMedianClusters;

for  i = 1:numel(SmallerMedianClusters(1,:))
    if SmallerMedianClusters(1,i) == 1
        figure(1); plot(mean(clWaveforms{i,2},2))
        title('clusters reducing median')
        hold on
    else disp('false')
       figure(2);plot(mean(clWaveforms{i,2},2))
       title('clusters increasing median')
       hold on
        
    end
    
end
msgbox('operation completed')

for i=1:size(clWaveforms) %indice con el rango
    p=cdfplot(test6(:,i)) %statement ( que tiene que hacer), pero se sobreescribe la variable en cada loop y
    %eso no lo queremos
    r(:,i)=p %creamos una variable contenedora de cada loop.
    hold on
end