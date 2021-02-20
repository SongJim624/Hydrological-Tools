function E = Budyko_Fu(AI, m)
    E = 1 + AI - (1 + AI .^ m) ^ (1 ./ m);
end