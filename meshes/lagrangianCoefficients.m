function lagrangCoeff = lagrangianCoefficients(indx, P, T, orderOfOccuracy)

nBasis = 0.5 * orderOfOccuracy * (orderOfOccuracy + 1);
coords = barycentr(orderOfOccuracy);

x = zeros(nBasis, 1);
y = zeros(nBasis, 1);

for i = 1 : nBasis
    for j = 1 : 3
        x(i) = x(i) + P(T(indx, j), 1) * coords(i,j);
        y(i) = y(i) + P(T(indx, j), 2) * coords(i,j);
    end
end

Vand = zeros(nBasis, nBasis);

k = 1;
for order = 0 : orderOfOccuracy-1
    for i = 0 : order
        for j = 1 : nBasis
            Vand(j,k) = Vand(j,k) + x(j) ^ (order-i) * y(j) ^ i;
        end
        k = k + 1;
    end
end

lagrangCoeff = inv(Vand);
