path = 'NACA0012_1_25deg/';
ax = [-0.5 1.25 -0.7 1.25];

% path = 'NACA0012_0deg/';
% ax = [-0.5 2 -1.5 1.5];
% ax = [-0.1 0.2 -0.15 0.15];

% path = 'GAMM-refined/';
% path = 'GAMM-doubleRefined/';
% path = 'channel';   
% ax = [0 3 0 1];

P = importdata(strcat(path, '/points.txt'));
T = importdata(strcat(path, '/triangles.txt'));
T = T + 1;

figure(7);
triplot(T, P(:,1), P(:,2), 'color', [0 0 0]);
set(gca, 'fontsize', 13);
% grid on
axis equal;   
axis(ax);

Solution = importdata(strcat(path, '/W.txt'));
W = solutionAtNodes(Solution, P, T);

% Solution = importdata(strcat(path, '/W-centres.txt'));
% W = calculateSolutionAtNodes(Solution', T, length(P));

rho = W(1,:)';
u =   W(2,:)' ./ rho;
v =   W(3,:)' ./ rho;
E =   W(4,:)';

kapa = 1.4;
p = (kapa-1) * (E - 1/2*rho.*(u.^2 + v.^2));
a = sqrt(kapa * p ./ rho);
mach = sqrt(u.^2 + v.^2) ./ a;

% colour of mesh
grey = [0.85 0.85 0.85];

% % % Mach contour plot % % %
figure(8);
% triplot(T, P(:,1), P(:,2), 'Color', grey);
% hold on;
[hC hC] = tricontf(P(:,1), P(:,2), T, mach, 30);
% set(hC,'LineStyle','none');
axis equal
axis(ax);
h = colorbar;
set(gca, 'fontsize', 13);
set(h, 'fontsize', 13);
% % % print(strcat(path,'mach3'), '-dpng', '-r800');

% grid on;

% artViscosity = importdata(strcat(path, 'W-artificialViscosity.txt'));
% hold on;
% for i = 1 : nTriangles
%     if (artViscosity(i) > 1e-4)
%         centre = 1/3 * (P(T(i,1), :) + P(T(i,2), :) + P(T(i,3), :));
%         scatter3(centre(1), centre(2), 1);
%     end
% end

% % % pressure contour plot % % %
figure(9);
% triplot(T, P(:,1), P(:,2), 'Color', grey);
% hold on;
[hC hC] = tricontf(P(:,1), P(:,2), T, p, 30);
% set(hC,'LineStyle','none');
axis equal
axis(ax);
h = colorbar;
set(gca, 'fontsize', 13);
set(h, 'fontsize', 13);
% % % print(strcat(path,'detailTlak2'), '-dpng', '-r800');

