function a = gen_general_partial_sequence(N, P, alpha, b)
a = zeros(N, 1);
for n = 0:N-1
    r = mod(n, P);
    s = (n-r)/P;
    phi = 2*pi/(N)*P*alpha*(n*(s+b)-P*s*(s+1)/2);
    a(n+1) = exp(1j*phi);
end
end
