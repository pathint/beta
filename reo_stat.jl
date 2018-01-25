function CountREO(A, B)
    if size(A) == size(B) || throw(DimensionMismatch("A does not match with B in size"))
        sum(A .< B)
     end
end

function MatrixREO(M)
    m, n = size(M) # m = sample size, n = genome size
    reos = zeros(Int64, m+1) 
    for i = 1:n-1
        for j = i+1:n
            reos[CountREO(M[:, i], M[:, j]) + 1] += 1
        end
    end
    return reos
end


# Example Usage
for i = 1:9
    # Expression profiles: m x n matrix (m, sample size; n genome size)
    data = readdlm(string("GSE92332_profile_", i, ".dat"), '\t')
    @time result = MatrixREO(data)
    writedlm(string("GSE92332_result_", i,".txt"), result, '\t')
end

