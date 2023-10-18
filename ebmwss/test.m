addpath('../gsc')
N = 64;
alpha = 1/5;
u = 0;
% bmwss
w1 = gen_ebmwss_sequence(N, alpha, u);
% gsc
P = 8;
b = 1/2;
w2 = gen_gsc_sequence(N, alpha, u, P, b);

W = {w1, w2};
thetas = linspace(0, 2*pi, 360);
A = zeros(length(thetas), length(W));
for ii = 1:length(W)
    w = W{ii};
    Nt = length(w);
    for i = 1:length(thetas)
        theta = thetas(i);
        arrayFactor = aTheta(Nt, cos(theta));
        A(i, ii) = w*arrayFactor';
        A(i, ii) = abs(A(i, ii))^2;
    end
end

close all
figure
polarplot(thetas, A, 'LineWidth', 1);
legend('EBMWSS', 'GSC')



function a = aTheta(N, shift)
a = exp(1j*pi*shift*(0:N-1));
end