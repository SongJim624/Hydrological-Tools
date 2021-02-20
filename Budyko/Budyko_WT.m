function E = Budyko_WT(AI, p)
    E = (1 + AI - sqrt((1 + AI) .^ p - 4 .* p .* (2- p) .* AI)) ./ (2 .* p .* (2 - p));
end