R = 1.1;
r = 0.9;

N = 600;
th = linspace(0, 2*pi, N);
circ = [cos(th); sin(th)];


c1 = [242, 101, 19]/255;
c2 = [255, 225, 79]/255;

h = 0.5;
% b = -0.2;
b = 0;
w = 0.3;
sp = 1;
% sp = 1.5;
arrow = [0, w, 0, -w, 0;
    h, b, 0, b, h];

arrowsh = [0; sp] + arrow;

% Nray = 15;
Nray = 16;
G = givens(cos(2*pi/Nray), sin(2*pi/Nray));


figure(1)
clf
hold on
patch(R*cos(th), R*sin(th), c1, 'Edgecolor', 'none')
patch(r*cos(th), r*sin(th), c2, 'Edgecolor', 'none')

for i = 1:Nray

    Gc = G^(i-1);

    arrowc = Gc*arrowsh;

    patch(arrowc(1, :), arrowc(2, :), c1, 'Edgecolor', 'none');

end

axis square
axis off