function [] = DrawROI(x, y, Color, Thickness, Linestyle)
%   draw a roi given by coordinates x,y 
if nargin < 3
    Color=[0.9, 0, 0];
    Thickness = 0.5;
    Linestyle = '-';
end
if nargin < 4
    Thickness = 0.5;
    Linestyle = '-';
end

if nargin < 5
    Linestyle = '-';
end

line(x, y, 'color', Color, 'LineWidth', Thickness, 'LineStyle', Linestyle);

end
