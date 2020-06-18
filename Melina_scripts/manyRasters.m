function [R lastlevel]=manyRasters(Spikes,col,h,ds,shift)

R={};
for i=1:numel(Spikes)
    R{i}=easyraster(Spikes{i},i+shift,col,h,ds);hold on
end
lastlevel=i+shift
ylim([0 lastlevel+1])



function r=easyraster(spikes,level,col,h,ds)

%this function produces a raster plot of vertical lines for multiple
%spike trains, spikes, a structure;
%Each spike train is plotted at increasing vertical coordinates, with height h,
%in color specified by string col, e.g. 'k','r'
%r is the final handle produced, this can be used to turn the baseline off
%(deals with bizarre behavior of Matlab 2010 stem plot).

r=[];
if ~isempty(spikes)
    if (size(spikes,2)<size(spikes,1)), spikes=spikes'; end %flip to row
    xs=[repmat(spikes,2,1);nan(size(spikes))];  %x coordinates for rasters
    xs=xs(:)/ds;
    ys=repmat([0;h;nan],numel(spikes),1);
    r=plot(xs,ys+level);
    set(r,'color', col, 'marker','none','linewidth',1);
end


