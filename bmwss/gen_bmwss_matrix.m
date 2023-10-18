function W = gen_bmwss_matrix(N1, N2, alpha, beta, u, v)
[w1, NA1] = gen_bmwss_sequence(N1, alpha, u);
[w2, NA2] = gen_bmwss_sequence(N2, beta, v);
W = w1.' * w2;
end