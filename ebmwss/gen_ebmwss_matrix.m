function W = gen_ebmwss_matrix(N1, N2, alpha, beta, u, v)
w1 = gen_ebmwss_sequence(N1, alpha, u);
w2 = gen_ebmwss_sequence(N2, beta, v);
W = w1.' * w2;
end