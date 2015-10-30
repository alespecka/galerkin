function NodeSolution = calculateSolutionAtNodes(Solution, T, nNodes)

[dim nTriangles nTimeSteps] = size(Solution);
NodeSolution = zeros(dim, nNodes, nTimeSteps);

for k = 1 : nTimeSteps
    NValues = zeros(nNodes, 1);
    for i = 1 : nTriangles
        for j = 1 : 3
            NodeSolution(:, T(i,j), k) = NodeSolution(:, T(i,j), k) + Solution(:, i, k);
            NValues(T(i,j)) = NValues(T(i,j)) + 1;
        end
    end
    
    for i = 1 : nNodes
        NodeSolution(:, i, k) = NodeSolution(:, i, k) / NValues(i);
    end
end

end
