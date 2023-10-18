function w = gen_gsc_sequence(N, alpha, u, P, b)
w = gen_general_partial_sequence(N, P, alpha, b);
shift = u-alpha;
w = w.' .* exp(1j*pi*(0:N-1)*shift);
w = w/norm(w);
end