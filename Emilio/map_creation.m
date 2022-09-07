CC_key = '';
for cc = 1:length(consCondNames)-1
    CC_key = [CC_key, sprintf('%s-', consCondNames{cc})]; %#ok<AGROW> 
end
CC_key = [CC_key, sprintf('%s', consCondNames{cc+1})];
CC_map = containers.Map(CC_key, {wruIdx});

C_key = condNames{chCond};
C_map = containers.Map(C_key, CC_map);

SW_key = sprintf('SW:%.2f-%.2f', spontaneousWindow*1e3);
SW_map = containers.Map(SW_key, C_map);

RW_key = sprintf('RW:%.2f-%.2f', responseWindow*1e3);
resp_map = containers.Map(RW_key, SW_map);