function b = barycentr(n)

x = linspace(0,1,n);
b = zeros(n*(n+1)/2, 3);
s = 1;
for i = 1:n
    for j = 1:i
        b(s,1) = 1-x(i);
        b(s,2) = x(j);
        b(s,3) = x(i) - x(j);
        s = s + 1;
    end
end
