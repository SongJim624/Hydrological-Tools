function E = Budyko_Old(AI)
    E = AI .* tanh(1 ./ AI);
end