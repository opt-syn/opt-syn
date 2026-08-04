function [fig] = theme_export(fname, fig, ext)
%THEME_EXPORT export a figure using light and dark themes

if nargin < 3
    ext = 'png';
end

if nargin < 2
    fig = gcf;
end

theme(fig, 'light');
f_light = [fname, '_light.', ext];
exportgraphics(fig, f_light);

theme(fig, 'dark');
f_dark = [fname, '_dark.', ext];
exportgraphics(fig, f_dark);

end