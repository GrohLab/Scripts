for x = 1:size(MatSignStep,2)
    for y = x+1:size(MatSignStep,1)
        fprintf(1, 'x:%d, y:%d\n', x, y)
    end
end