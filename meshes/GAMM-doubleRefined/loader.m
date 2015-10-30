load geometrie;

P = zeros(length(geometrie{1}), 2);
P(:,1) = geometrie{1};
P(:,2) = geometrie{2};
T = geometrie{3};
adjM = geometrie{4};
T = T(:,1:3);
adjM = adjM(:,1:3);

dlmwrite('points.txt', P, 'delimiter', ' '); %dlmwrite('myFile.txt',M,'delimiter','\t','precision',3)
dlmwrite('triangles.txt', T, 'delimiter', ' ');
dlmwrite('neighbours.txt', adjM, 'delimiter', ' ');
