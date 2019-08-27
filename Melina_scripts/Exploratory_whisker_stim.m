%% looking at whisker triggers practice script


figure
plot(t_supply_voltage,supply_voltage_data)

%%

close all
figure
%look at each channel
for i=1:16
    plot(board_dig_in_data(i,:))
    title(i)
    pause
end


%% plot puff and starts of puff
close all
puff=board_dig_in_data(4,:);

t_puff=t_dig(find(puff));


figure
plot(t_dig,puff);hold on
plot(t_puff,ones(size(t_puff)),'r.')
xlabel('seconds')

min_dt=5
dt_puff=[min_dt diff(t_puff)];
puff_starts=find(dt_puff>=min_dt);
plot(t_puff(puff_starts),ones(size(puff_starts)),'g*')


save WhiskerDataRough puff t_dig
%%
figure
histogram(log10(dt_puff*1000),500)

