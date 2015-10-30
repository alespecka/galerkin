function combinedW = linearCombination(w, x, y, centreX, centreY, powers, nBasis)

    combinedW = zeros(4, 1);
    for j = 1 : nBasis
        basisFun = (x - centreX) ^ powers(j,1) * (y - centreY) ^ powers(j,2);
        for d = 1 : 4
            combinedW(d) = combinedW(d) + w((d-1)*nBasis+j) * basisFun;
        end
    end
    
end
