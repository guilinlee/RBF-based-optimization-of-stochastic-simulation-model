function [h]=plot3d_errorbars(x, y, z, u, l, t)
%checked#4

h=plot3(x, y, z, '.k');

set(h, 'MarkerSize', 20);
hold on

for i=1:length(x)
        xV = [x(i); x(i)];
        yV = [y(i); y(i)];
        zV = [l(i), u(i)];
        h=plot3(xV, yV, zV, '-k');
        set(h, 'LineWidth', 2);
end
if(isempty(t)~=1)
if(size(t,2)==1)
    h=plot3(x, y, t, '*');
else
    h=plot3(t(:,1), t(:,2), t(:,3), '*');
end
end

set(h, 'MarkerSize', 15);
box on, grid on