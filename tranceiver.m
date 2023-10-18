function errorNum = tranceiver(bpsk, w, snr, alpha, beta, u0, v0, chan)
u1 = gen_rand_pos(u0, alpha);
v1 = gen_rand_pos(v0, beta);
[N1, N2] = size(w{1});
H = exp(-1j*pi*u1*(0:N1-1)).' * exp(-1j*pi*(v1*(0:N2-1)));
if isequal(chan, 'rayleigh')
    H = 1/sqrt(2)*(randn+1j*randn) * H;
end
N = length(bpsk);
bits = (bpsk>0);
noise = 1/sqrt(snr) * 1/sqrt(2)*(randn(1, N)+1j*randn(1, N));

if length(w) == 1
    y = abs(w{1}(:)'*H(:))* bpsk + noise;
else
    y = sqrt(0.5*(abs(w{1}(:)'*H(:))^2 + abs(w{2}(:)'*H(:))^2)) * bpsk + noise;
end
errorNum = sum(abs(bits-(real(y)>0)));
end

function u1 = gen_rand_pos(u0, alpha)
u1 = (2*rand-1)*alpha/1 + u0;
end

