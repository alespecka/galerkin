path = 'channel/';
P = importdata(strcat(path, 'points.txt'));
T = importdata(strcat(path, 'triangles.txt'));
W = importdata(strcat(path, 'nodeSolution.txt'));
T = T + 1;


figure(5);
hold on
for i = 1 : size(W, 1);
    rho = [W(i,1) W(i,5) W(i,9)];
    u = [W(i,2) W(i,6) W(i,10)] ./ rho;
    v = [W(i,3) W(i,7) W(i,11)] ./ rho;
    E = [W(i,4) W(i,8) W(i,12)];

    kapa = 1.4;
    p = (kapa-1) * (E - 1/2*rho.*(u.^2 + v.^2));
    a = sqrt(kapa * p ./ rho);
    mach = sqrt(u.^2 + v.^2) ./ a;

    trisurf([1 2 3], P(T(i,:),1)', P(T(i,:),2)', mach);
end

view(0,0);
%axis([0 3 0 1 0 1]);
    