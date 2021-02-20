function E = Budyko_Zhang(AI, w)
    E = (1 + w .* AI) ./ (1 + w .* AI + 1 ./ AI);
end