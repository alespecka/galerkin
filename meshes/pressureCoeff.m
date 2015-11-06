% function pressureCoeff(path, machInf)

machInf = 0.5;
path = 'NACA0012_0deg/';

P = importdata(strcat(path, '/points.txt'));
T = importdata(strcat(path, '/triangles.txt'));
N = importdata(strcat(path, '/neighbours.txt'));
T = T + 1;
N = N + 1;

kapa = 1.40;
nTriangles = size(T, 1);
nFaces = 3;


% Solution = importdata(strcat(path, '/W.txt'));
W = importdata(strcat(path, '/W-centres.txt'));

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

% p = (kapa-1) * (W(:,4) - 1/2*(W(:,2).^2 + W(:,3).^2) ./ W(:,1));
% 
% pInf = (1 + (kapa-1)/2 * machInf^2) ^ ((kapa-1)/kapa);
% rhoInf = (1 + (kapa-1)/2 * machInf^2) ^ (1/(kapa-1));
% uInf = machInf * sqrt( (kapa * pInf) / rhoInf);

nVals = 0;
for i = 1 : nTriangles    
    for j = 1 : nFaces
        if N(i,j) == 0 && P(T(i,j),2) >= 0
            nVals = nVals + 1;
        end
    end
end

cp = zeros(nVals, 1);
x = zeros(nVals, 1);
pos = 0;
for i = 1 : nTriangles
    for j = 1 : nFaces
        if N(i,j) == 0 && P(T(i,j),2) >= 0
            pos = pos + 1;
            cp(pos) = 2 * (p(i) - pInf) / (rhoInf * uInf^2);
            
            x1 = P(T(i,j), 1);
            x2 = P(T(i,mod(j,3)+1), 1);
            x(pos) = (x1 + x2) / 2;
        end
    end
end

[x, I] = sort(x);
cp = cp(I);

figure(1);
plot(x, cp, '-', 'Color', [0 0 0], 'linewidth', 2)
set(gca,'Ydir','reverse');
set(gca, 'fontsize', 15);
xlabel('x/c');
ylabel('c_p');
% L = get(gca,'YLim');
axis([0 1 -0.5 1])
% set(gca,'YTick', linspace(L(1),L(2), 6))
grid on

% print(strcat(path,'tlakKoef'), '-dpng', '-r800');