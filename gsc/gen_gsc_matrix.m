function W = gen_gsc_matrix(N1, N2, alpha, beta, u, v)
P1 = 1;
P2 = 1;
b = 1/2;
w1 = gen_gsc_sequence(N1, alpha, u, P1, b);
w2 = gen_gsc_sequence(N2, beta, v, P2, b);
W = w1.' * w2;
end