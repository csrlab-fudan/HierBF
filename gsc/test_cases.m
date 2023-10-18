%% Milewski sequence (same length)
M = 3;
k = 2;
N = M^(2*k+1);
P = M^k;
alpha = 1;
b = 1;
a = gen_general_partial_sequence(N, P, alpha, b);
close all
figure
stem(abs(fft(a)).^2)

%% perfect sequence
P = 4;
Q = 5;
N = P^2*Q;
alpha = 1;
if mod(P+Q, 2)==1
    b = 1/2;
else
    b = 1;
end
a = gen_general_partial_sequence(N, P, alpha, b);
omega0 = exp(1j*2*pi/P);
omega1 = exp(1j*2*pi/Q);
omega2 = exp(1j*2*pi/(P*Q));
omega3 = exp(1j*2*pi/(N));
f = zeros(1, N);
for i = 0:N-1
    v = mod(i, P);
    u = (i-v)/P;
    if mod(P, 2)==1
%         for n = 0:N-1
%             f(i+1) = f(i+1) + omega3^(-n*i)*a(n+1);
%         end

%         for s = 0:P*Q-1
%             for r = 0:P-1
%                 f(i+1) = f(i+1) + omega3^(-(s*P+r)*((u*P+v)-P*(s+b))-0.5*P^2*s*(s+1));
%             end
%         end

%         for s = 0:P*Q-1
%             for r = 0:P-1
%                 f(i+1) = f(i+1) + omega1^(s*((s-1)/2+b-u))*omega2^(s*(r-v)+r*(b-u))*omega3^(-r*v);
%             end
%         end

%         for r = 0:P-1
%             for w = 0:Q-1
%                 for t = 0:P-1
%                     f(i+1) = f(i+1) + (-1)^(t*(Q+2*b-1))*omega0^(t*(r-v))*omega1^(0.5*w^2+(b-0.5-u)*w)*...
%                         omega2^(w*(r-v)+r*(b-u)) * omega3^(-r*v);
%                 end
%             end
%         end

%         f(i+1) = P*omega2^(v*(b-u))*omega3^(-v^2)*sum(omega1.^(0.5*(0:Q-1).^2+(b-0.5)*(0:Q-1)));
        f(i+1) = P*omega2^(v*(b-u))*omega3^(-v^2)*sqrt(Q)*exp(1j*pi*(1/4-(2*b-2*u-1)^2/(4*Q)));
    elseif v<P/2
        f(i+1) = P*omega2^((v+P/2)*(b-u))*omega3^(-v*(v+P/2))*sqrt(Q)*exp(1j*pi*(1/4-(b-u)^2/Q));
    else
        f(i+1) = P*omega2^((v-P/2)*(b-u))*omega3^(-v*(v-P/2))*sqrt(Q)*exp(1j*pi*(1/4-(b-u-1)^2/Q));
    end
end
close all
figure
stairs(angle(f))
hold on
stairs(angle(fft(a)))
legend('exact', 'fft')

figure
% stem(abs(f))
% hold on
stem(abs(fft(a, 2*N)))
% legend('exact', 'fft')

%% perfect sequence
P = 5;
Q = 2;
N = P^2*Q;
alpha = 1;
if mod(P+Q, 2)==1
    b = 1/2;
else
    b = 1;
end
a = gen_general_partial_sequence(N, P, alpha, b);

close all
fftLen = N;
stem(abs(fft(a, fftLen)))
