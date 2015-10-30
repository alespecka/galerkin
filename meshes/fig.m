% path = 'NACA0012b/';
% ax = [-0.4 1.1 -0.5 1];

% path = 'channel/';
ax = [0 3 0 1];

P = importdata(strcat(path, '/points.txt'));
T = importdata(strcat(path, '/triangles.txt'));
T = T + 1;
    
W = calculateSolutionAtNodes(Solution', T, length(P));

rho = W(1,:)';
u = W(2,:)' ./ rho;
v = W(3,:)' ./ rho;
E = W(4,:)';

kapa = 1.4;
p = (kapa-1) * (E - 1/2*rho.*(u.^2 + v.^2));
a = sqrt(kapa * p ./ rho);
mach = sqrt(u.^2 + v.^2) ./ a;

figure(7);    
triplot(T, P(:,1), P(:,2));
axis([-2.5, 3.5, -3, 3]);
set(gca, 'fontsize', 13);
grid on
axis equal;   
axis(ax);
% print('nesmysl', '-dpng', '-r800');

% Solution = importdata(strcat(path,'/centres.txt'));
Solution = importdata(strcat(path, '/W-centres.txt'));

figure(8);
tricontf(P(:,1), P(:,2), T, mach, 30);
axis equal
axis(ax);
h = colorbar;
set(gca, 'fontsize', 13);
set(h, 'fontsize', 13);

figure(9);
tricontf(P(:,1), P(:,2), T, p, 30);
axis equal
axis(ax);
h = colorbar;
set(gca, 'fontsize', 13);
set(h, 'fontsize', 13);

% figure(4);
% trisurf(T, P(:,1), P(:,2), mach);
% view(0,0);
% axis([0 3 0 1     0 1]);

% figure(4);
% trisurf(T, P(:,1), P(:,2), mach);
% view(0,90);

% hold on;
% 
% X = importdata('X.txt');
% Y = importdata('Y.txt');
% 
% for i = 1 : length(X)
%     scatter3(X(i), Y(i), 1.5);
% end
