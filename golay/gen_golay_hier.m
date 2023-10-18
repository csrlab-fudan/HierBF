function [W1, W2] = gen_golay_hier(N1, N2, Ns1, Ns2, u, v)
[omni1, omni2] = gen_quaternary_acm_pair(N1/Ns1, N2/Ns2);
omni1 = omni1/norm(omni1(:));
omni2 = omni2/norm(omni2(:));
a1 = 1/sqrt(Ns1)*exp(1j*pi*(0:Ns1-1)*u);
a2 = 1/sqrt(Ns2)*exp(1j*pi*(0:Ns2-1)*v);
direct = a1.'*a2;
W1 = kron(omni1, direct);
W2 = kron(omni2, direct);
end
