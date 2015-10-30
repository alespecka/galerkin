function combinedW = linearCombinationLagrange(w, x, y, lagrangeCoeff, orderOfOccuracy)

    combinedW = zeros(4, 1);
    nBasis = 0.5 * orderOfOccuracy * (orderOfOccuracy + 1);
    
    for i = 1 : nBasis
        basisFun = 0;
        k = 1;
        for order = 0 : orderOfOccuracy-1
            for j = 0 : order
                basisFun = basisFun + lagrangeCoeff(k,i) * x ^ (order-j) * y ^ j;
                k = k + 1;
            end
        end

        for d = 1 : 4
            combinedW(d) = combinedW(d) + w((d-1)*nBasis+i) * basisFun;
        end
    end
    
end