Nm = numel(mice);
Ns = arrayfun(@(m) numel( m.Sessions ), mice );
aiPop = zeros( sum( Ns ), 4 );
aiPbp = zeros( sum( Ns ), 8, 2 );
aiPbp2 = zeros( 8, 2, sum( Ns ) );
cr = 1;
for cm = 1:Nm
    for cs = 1:Ns(cm)
        dt = mice(cm).Sessions(cs).DataTable;
        aiPop(cr,:) = [cm, cs, dt.AmplitudeIndex];
        aiPbp(cr,:,:) = reshape( dt.AI_pbp{:}, 1, [], 2 );
        aiPbp2(:,:,cr) = dt.AI_pbp{:};
        cr = cr + 1;
    end
end