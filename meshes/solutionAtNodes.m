function NodeSolution = solutionAtNodes(W, P, T)

nEqns = 4;
nFaces = 3;

nNodes = size(P, 1);
nTriangles  = size(T, 1);
nBasis = size(W,2) / nEqns;

orderOfOccuracy = 0;
n = nBasis;
while n > 0
    n = n - orderOfOccuracy - 1;
    orderOfOccuracy = orderOfOccuracy + 1;
end

NodeSolution = zeros(nEqns, nNodes);
NValues = zeros(nNodes, 1);
for i = 1 : nTriangles
    lagrangeCoeffs = lagrangianCoefficients(i, P, T, orderOfOccuracy);
    
    for j = 1 : nFaces
        x = P(T(i,j), 1);
        y = P(T(i,j), 2);
        combinedW = linearCombinationLagrange(W(i,:), x, y, lagrangeCoeffs, orderOfOccuracy);
        NodeSolution(:, T(i,j)) = NodeSolution(:, T(i,j)) + combinedW;
        NValues(T(i,j)) = NValues(T(i,j)) + 1;
    end
end

for i = 1 : nNodes
    NodeSolution(:, i) = NodeSolution(:, i) / NValues(i);
end

