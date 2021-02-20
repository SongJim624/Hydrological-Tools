function E = Budyko(AI)
E = sqrt(AI .* tanh(1 ./ AI) .* (1 - exp(-AI)));
end