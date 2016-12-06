function plotCircle(x0,r,colorstring)

for n = 1:40
    phi = (n-1)*2*pi / (40-1);
    u = [cos(phi) sin(phi)]';
    x(:,n) = x0 + r * u;
end

hold on;
plot(x(1,:) + i*x(2,:), colorstring, 'lineWidth', 2);

return

