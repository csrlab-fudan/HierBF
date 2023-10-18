function [w, NA] = gen_bmwss_sequence(N, alpha, u)
NA = [];
for k = 0:N-1
    if mod(sqrt((N-k)*alpha), 1)==0 && mod(N-k, sqrt((N-k)*alpha))==0
        NA = N-k;
        break
    end
end
if isempty(NA)
    error('Infeasible length!')
end

M = sqrt(NA*alpha);
w = gen_phase_shift(N, M, NA, u);
end

function w = gen_phase_shift(N, M, NA, u)
Ns = NA/M;
if mod(Ns, 2)==0
    phaseShift = exp(-1j*(1:M)*(Ns-1)/Ns*pi);
else
    phaseShift = exp(1j*(1:M)*pi/Ns);
end

w = zeros(1,N);
for i = 1:M
    w(1+(i-1)*Ns:i*Ns) = phaseShift(i)*aTheta(Ns, -1+(2*i-1)/Ns);
end
shift = 1-M^2/NA+u;
w = w.*aTheta(N, shift);
w = w/sqrt(N);
end

function a = aTheta(N, shift)
a = exp(1j*pi*shift*(0:N-1));
end