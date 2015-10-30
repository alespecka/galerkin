function [areas indiameters centres lens normals] = triangleProps(P, T)

nTriangles = size(T, 1);
lens = zeros(nTriangles, 3);
normals = zeros(nTriangles, 6);
areas = zeros(nTriangles, 1);
indiameters = zeros(nTriangles, 1);
centres = zeros(nTriangles, 2);

for i = 1 : nTriangles
    points = T(i, [1 2 3]);

    for j = 1 : 3
        a = points(j);
        b = points(mod(j,3) + 1); 
        n = [P(b,2)-P(a,2) P(a,1)-P(b,1)]';
        lens(i,j) = norm(n);
        normals(i, [2*j-1 2*j]) = n / lens(i,j);
    end

    M = [P(T(i,1),1) P(T(i,1),2);
         P(T(i,2),1) P(T(i,2),2);
         P(T(i,3),1) P(T(i,3),2)];

    centres(i,:) = mean(M);
    
    areas(i) = 1/4 * sqrt((lens(i,1)+lens(i,2)+lens(i,3)) * ...
        (-lens(i,1)+lens(i,2)+lens(i,3)) * (lens(i,1)-lens(i,2)+lens(i,3)) * ...
        (lens(i,1)+lens(i,2)-lens(i,3)));

    indiameters(i) = 4 * areas(i) / (lens(i,1)+lens(i,2)+lens(i,3));
end

end
