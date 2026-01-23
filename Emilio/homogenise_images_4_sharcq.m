
fullpath = @(x) fullfile(x.folder, x.name);
get_img_size = @(x) [x.getTag("ImageWidth"), x.getTag("ImageLength")];
fnOpts = {'UniformOutput', false};
ideal_prop = 57/40;


img_dir = "Z:\Leonie\AC\M65\M65 ThSC";
original_resolution = 2.302;
desired_resolution = 10;
homogenised_img_folder = fullfile(img_dir, 'Homogenised and downsampled');
if ~exist(homogenised_img_folder,"dir")
    mkdir(homogenised_img_folder)
end
scaling_factor = original_resolution/desired_resolution;

img_paths = dir(fullfile(img_dir, "*Merged_overlay.tif"));
roi_paths = dir(fullfile(img_dir, "*Merged_overlay ROIs.csv"));
t_objs = arrayfun(@(x) Tiff(fullpath(x), 'r'), img_paths);
img_size = arrayfun(@(x) get_img_size(x), t_objs, fnOpts{:});
arrayfun(@(x) x.close, t_objs);
img_size = cat(1, img_size{:});
img_prop = img_size(:,1)./img_size(:,2);
final_sz = ceil(max(img_size(:,1)) * [1, (1/ideal_prop)]);
for ci = 1:numel(img_paths)
    shift_px = [floor((final_sz([2,1]) - img_size(ci,[2,1]))/2)];
    img = imread(fullpath(img_paths(ci)));
    img2 = padarray(img, [shift_px, 0], 0, 'both');
    img2 = imresize(img2, scaling_factor);
    roi_table = readtable(fullpath(roi_paths(ci)), ...
        'Delimiter',',','VariableNamingRule','preserve');
    roi_table{:,{'X','Y'}} = round((roi_table{:,{'X','Y'}} + ...
        shift_px([2,1]))*scaling_factor);
    imwrite(img2,fullfile(homogenised_img_folder, ...
        img_paths(ci).name),'tif', 'WriteMode', 'overwrite')
    writetable(roi_table,fullfile(homogenised_img_folder, ...
        roi_paths(ci).name),'WriteVariableNames',true)
end