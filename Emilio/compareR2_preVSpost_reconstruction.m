mice_results = "Z:\Emilio\SuperiorColliculusExperiments\Roller\PoolFigures\MC-iegRNs";
load( fullfile( mice_results, 'iRNs\MCiRNs_reconstruction_sm.mat') )

Nm = numel( mice ); % Numer of mice
Nspm = arrayfun(@(x) numel( x.Sessions ), mice ); % Number of sessions per mouse
Nexp = sum( Nspm );
r2_res_c = zeros( 3, 8, Nexp );
r2_res_l = zeros( 3, 8, Nexp );
ce = 1;
mouseID = zeros( Nexp, 1 );
sessID = zeros( Nexp, 1 );
for cm = 1:Nm
    for cs = 1:Nspm(cm)
        dt = mice(cm).Sessions(cs).DataTable;
        r2 = dt.R_2_p_L;
        r2_res_c(:,:,ce) = r2{1};
        r2_res_l(:,:,ce) = r2{2};
        mouseID(ce) = cm;
        sessID(ce) = cs;
        ce = ce + 1;
    end
end