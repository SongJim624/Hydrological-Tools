function E = Budyko_CY(AI, n)
    E = 1 ./ ((1 ./ AI) ^ n + 1) .^ n;
end