%subdifferential of an indicator function

FS = 14;
FST = 20;

h = 1;
b = 5;
% b = 10;

Gra_F = blkdiag(h, b) * ...
    [-1, -1, 1, 1;
    -1, 0, 0, 1];

J = [0, 1; 1, 0];
Gra_Finv = J * Gra_F;

%% plot F and Finv
figure(4)
clf;
tiledlayout(1, 2)
nexttile
hold on
plot(Gra_F(1, :), Gra_F(2, :), 'k', 'linewidth', 4);
plot([-b, b], [0,0], ':k')
plot([0,0], [-b, b], ':k')
xlabel('$z$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$w$', 'interpreter', 'latex', 'fontsize', FS)
title('$w \in F(z)$', 'interpreter', 'latex', 'fontsize', FST)
xlim([-b, b])
ylim([-b, b])
axis square

nexttile
hold on
plot(Gra_Finv(1, :), Gra_Finv(2, :), 'k', 'linewidth', 4);
plot([-b, b], [0,0], ':k')
plot([0,0], [-b, b], ':k')
xlabel('$w$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$z$', 'interpreter', 'latex', 'fontsize', FS)
title('$z \in F^{-1}(w)$', 'interpreter', 'latex', 'fontsize', FST)
xlim([-b, b])
ylim([-b, b])
axis square

% pause

%% plot resolvent and yosida
figure(4)
clf;
tiledlayout(1, 2)


BAD = false;
if BAD
    be = 1;
    bsign = '';
else
    be = -1;
    bsign = '-';
end

Elist = [0.1, 0.5, 1, 2, 10];
sname = {'\frac{1}{10}', '\frac{1}{2}', '1', '2', '10'};



Ne = length(Elist);    

nexttile
hold on
for i = 1:Ne
    E = be*Elist(i);
    Tw = [1, -E; %v = z - E w
      0, 1]; %w = w
    Gra_w = Tw * Gra_F;
    plot(Gra_w(1, :), Gra_w(2, :), 'linewidth', 4);
end
plot([-b, b], [0,0], ':k')
plot([0,0], [-b, b], ':k')
xlabel('$v$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$w$', 'interpreter', 'latex', 'fontsize', FS)
title('$w \in (F^{-1} - E)^{-1}(v)$', 'interpreter', 'latex', 'fontsize', FST)
xlim([-b, b])
ylim([-b, b])

axis square
lname = cell(Ne, 1);
for i = 1:Ne
    lname{i} = ['$E=', bsign, sname{i}, '$'];
end
legend(lname, 'location', 'southeast', 'interpreter', 'latex', 'fontsize', FS*0.8)


nexttile
hold on
for i = 1:Ne
    E = be*Elist(i);
    Tz = [1, -E; %v = z - E w
         1, 0];  %z = z
            
    Gra_z = Tz * Gra_F;
    

    plot(Gra_z(1, :), Gra_z(2, :),  'linewidth', 4);
end
plot([-b, b], [0,0], ':k')
plot([0,0], [-b, b], ':k')
xlabel('$v$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$z$', 'interpreter', 'latex', 'fontsize', FS)
title('$z \in (I - EF)^{-1}(v)$', 'interpreter', 'latex', 'fontsize', FST)
xlim([-b, b])
ylim([-b, b])
axis square

%% sweep for incremental positivity

N = 100;
w1 = linspace(-b, 0, N);
z1 = -h*ones(1, N);

z2 = linspace(-h, h, N);
w2 = zeros(1, N);


w3 = linspace(0, b, N);
z3 = h*ones(1, N);

w = [w1, w2, w3];
z = [z1, z2, z3];

%incremental differences
dz = bsxfun(@minus,z,z');
dw = bsxfun(@minus,w,w');

% figure(6)
% clf
% scatter(dz, dw, 'k')