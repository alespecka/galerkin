path = '../NACA0012d/';
P = importdata(strcat(path, 'points.txt'));
T = importdata(strcat(path, 'triangles.txt'));
T = T + 1;
N = importdata(strcat(path, 'neighbours.txt'));
N = N + 1;
W = importdata(strcat(path, 'centres.txt'));
% W = importdata('../centres.txt');

rho = W(:,1);
u = W(:,2) ./ rho;
v = W(:,3) ./ rho;
E = W(:,4);

kapa = 1.4;
machInf = 0.5;
p = (kapa-1) * (E - 1/2*rho.*(u.^2 + v.^2));
a = sqrt(kapa * p ./ rho);
mach = sqrt(u.^2 + v.^2) ./ a;

pInf = 1 / (1 + (kapa-1)/2 * machInf^2) ^ (kapa/(kapa-1));
rhoInf = (1 + (kapa-1)/2 * machInf^2) ^ (1/(1-kapa));
uInf = machInf * sqrt(kapa*pInf/rhoInf);

% [I,J] = find(N == -3);

nTriangles = length(T);
I = [];
J = [];

[areas, indiameters, centres, lens, normals] = triangleProps(P, T);

for i = 1 : nTriangles;
    if abs(centres(i,1)) < 3 && abs(centres(i,2)) < 3
        for j = 1 : 3
            if N(i,j) == 0 || N(i,j) == -3
                I = [I,i];
                J = [J,j];
            end
        end
    end
end

px = zeros(length(I), 1);
py = zeros(length(I), 1);
cp = zeros(length(I), 1);
x12 = zeros(length(I), 1);

% [areas, indiameters, centres, lens, normals] = triangleProps(P, T);
% max = 0; min = 1000;
for i = 1 : length(I)
    x1 = P(T(I(i),J(i)), 1);
%     y1 = P(T(I(i),J(i)), 2);
    x2 = P(T(I(i),mod(J(i),3)+1), 1);
%     y2 = P(T(I(i),mod(J(i),3)+1), 2);
%     
    x12(i) = (x2 + x1) / 2;
%     y12 = (y2 + y1) / 2;
%     
%     if max < y12
%         max = y12;
%     end
%     
%     if min > y12
%         min = y12;
%     end
%     
%     x = x2 - x1;
%     y = y2 - y1;
%     len = sqrt(x^2 + y^2);
%     
%     nx = y / len;
%     ny = -x / len;
    
%     px(i) = p(I(i)) * nx * len;
%     py(i) = p(I(i)) * ny * len;

    cp(i) = 2 * (p(I(i)) - pInf) / (rhoInf * uInf^2);
    
    px(i) = p(I(i)) * normals(I(i),2*J(i)-1) * lens(I(i), J(i));
    py(i) = p(I(i)) * normals(I(i),2*J(i)) * lens(I(i), J(i));
end

% A = max - min;

cdx = 2 * sum(px) / (rhoInf * uInf^2);
cdy = 2 * sum(py) / (rhoInf * uInf^2);

figure(10);
scatter(x12, cp, '.');

figure(11);
scatter(x12, p(I), '.');

figure(12);
scatter(x12, mach(I), '.');

