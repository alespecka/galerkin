path = 'channel/';
% path = 'NACA0012d/';

P = importdata(strcat(path, '/points.txt'));
T = importdata(strcat(path, '/triangles.txt'));
W = importdata(strcat(path, '/W.txt'));
T = T + 1;

nTriangles  = size(W, 1);
nBasis = size(W,2) / 4;

orderOfOccuracy = 0;
n = nBasis;
while n > 0
    n = n - orderOfOccuracy - 1;
    orderOfOccuracy = orderOfOccuracy + 1;
end

n = 1;
powers = zeros(nBasis, 2);
for degree = 0 : orderOfOccuracy - 1
    for j = 0 : degree
        powers(n,1) = degree - j;
        powers(n,2) = j;
        n = n + 1;
    end
end

[areas, indiameters, centres, lens, normals] = triangleProps(P, T);

kapa = 1.4;

b = barycentr(8);
N = size(b,1);
% N = 4;
p = zeros(N, 1);
mach = zeros(N, 1);
x = zeros(N, 1);
y = zeros(N, 1);

x(j) = 0;
y(j) = 0;

for j = 1 : N
    for l = 1 : 3
        x(j) = x(j) + b(j,l) * P(T(1,l), 1);
        y(j) = y(j) + b(j,l) * P(T(1,l), 2);
    end
end
TRI = delaunay(x,y);

minMach = 100;
maxMach = 0;
figure(10);
hold on
for i = 1 : nTriangles
    lagrangeCoeffs = lagrangianCoefficients(i, P, T, orderOfOccuracy);
    for j = 1 : N
        x(j) = 0;
        y(j) = 0;
        for l = 1 : 3
            x(j) = x(j) + b(j,l) * P(T(i,l), 1);
            y(j) = y(j) + b(j,l) * P(T(i,l), 2);
        end
        
        combinedW = linearCombinationLagrange(W(i,:), x(j), y(j), lagrangeCoeffs, orderOfOccuracy);
%         combinedW = linearCombination(W(i,:), x(j), y(j), centres(i,1), centres(i,2), powers, nBasis);
        
        rho = combinedW(1);
        u = combinedW(2) ./ rho;
        v = combinedW(3) ./ rho;
        E = combinedW(4);

        p(j) = (kapa-1) * (E - 1/2*rho.*(u.^2 + v.^2));
        a = sqrt(kapa * p(j) ./ rho);
        mach(j) = sqrt(u.^2 + v.^2) ./ a;
        
        if mach(j) < minMach
            minMach = mach(j);
        end
        
        if mach(j) > maxMach
            maxMach = mach(j);
        end
    end
    
    trisurf(TRI, x, y, mach, 'edgecolor', 'none'); %'LineWidth',2
end
axis equal
% shading interp
% axis([0 3 0 1 0 2])
%view(0,0);
% shading interp

% contourLine = linspace(minMach, maxMach, 30)';
% 
% figure(10);
% hold on
% for i = 1 : nTriangles    
%     for j = 1 : N
%         x(j) = 0;
%         y(j) = 0;
%         for l = 1 : 3
%             x(j) = x(j) + b(j,l) * P(T(i,l), 1);
%             y(j) = y(j) + b(j,l) * P(T(i,l), 2);
%         end
%         
%         combinedW = linearCombination(W(i,:), x(j), y(j), centres(i,1), centres(i,2), powers, nBasis);
%         
%         rho = combinedW(1);
%         u = combinedW(2) ./ rho;
%         v = combinedW(3) ./ rho;
%         E = combinedW(4);
% 
%         p(j) = (kapa-1) * (E - 1/2*rho.*(u.^2 + v.^2));
%         a = sqrt(kapa * p(j) ./ rho);
%         mach(j) = sqrt(u.^2 + v.^2) ./ a;
%     end
%     TRI = delaunay(x,y);
%     tricontf(x, y, TRI, mach);
% end
% colorbar;

