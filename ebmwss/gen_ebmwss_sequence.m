function w = gen_ebmwss_sequence(N, alpha, u)
M = 0;
for i = 1:N
    if mod(N, i) == 0 && i^2>=N*alpha
        M = i;
        break
    end
end
assert(M~=0, 'Infeasible length!')
delta = alpha/M;
Ns = N/M;
if mod(Ns, 2)==0
    phaseShift = exp(-1j*pi*delta*(1:M).*(Ns*(1:M)-1));
else
    phaseShift = exp(-1j*pi*delta*Ns*(1:M).^2);
end

w = zeros(1,N);
for i = 1:M
    w(1+(i-1)*Ns:i*Ns) = phaseShift(i)*aTheta(Ns, u-alpha+(2*i-1)*delta);
end
w = w/sqrt(N);
% for i = 1:M-1
%     inter_point = u-alpha+2*i*delta;
%     display(abs(w*aTheta(N, inter_point)')^2)
% end
end

function a = aTheta(N, shift)
a = exp(1j*pi*shift*(0:N-1));
end