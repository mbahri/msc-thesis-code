function [] = plot_with_detail(O, xrange, yrange, zoom_factor)

[h, w, ~] = size(O);

if nargin < 4
    zoom_factor = 2;
end

% fig = figure;
image(O);
% imshow(O, []);
set(gca,'visible','off')
axes('position',[.65 .175 .25 .25])
box on
image(O(xrange, yrange, :));
% subimage(O(xrange, yrange, :));

axis tight
zoom(zoom_factor)
set(gca,'XtickLabel',[],'YtickLabel',[]);

fig.Units = 'pixels';
fig.Position = [0 0 w h];

end