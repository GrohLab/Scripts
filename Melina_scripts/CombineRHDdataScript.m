%% stitching files together

clear all
read_Intan_RHD2000_file  %to do is to remove gui and just do a loop to import all data in directory
display('saving...')
save(fullfile(path,filename(1:end-4)), '-v7.3')
display('done')



%% outline of strategy

%do a loop over all of the .mat files
    %from each file, read in only amplifier_data variable
    %concatenate amplifier_data     MAYBE WITH  >> cat
%save the output

%%