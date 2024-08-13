fnOpts = { 'UniformOutput', false };
Ni = arrayfun(@(m) arrayfun(@(s) size( s.Intensities, 1 ),  m.Sessions ), ...
    mice, fnOpts{:} );
Nr = sum( cat( 1, Ni{:} ) );
hab_table = zeros( Nr, 12 );
cr = 1;
for cm = 1:numel(mice)
    for cs = 1:numel(mice(cm).Sessions)
        Nism = size( mice(cm).Sessions(cs).Intensities, 1 );
        idxs = (0:Nism-1) + cr;
        hab_table(idxs, 1:2) = repmat( [cm,cs], Nism, 1 );
        hab_table(idxs, 3:4) = [mice(cm).Sessions(cs).Intensities, ...
            mice(cm).Sessions(cs).BehIndex];
        hab_table(idxs, 5:end) = mice(cm).Sessions(cs).MaxVals;
        cr = cr + Nism;
    end
end