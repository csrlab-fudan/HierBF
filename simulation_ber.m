N1 = 40;
addpath('gsc', 'golay', 'bmwss', 'ebmwss');
N2 = 36;
u = 0;
v = 0;
targetTheta = pi/2;
targetPhi = 0;
alpha = 1/4;
beta = 1/4;

[W1Golay, W2Golay] = gen_golay_hier(N1, N2, 1/alpha, 1/beta, u, v);
W1DFT = W1Golay/sqrt(N1*N2*alpha*beta);
W2DFT = W1Golay/sqrt(N1*N2*alpha*beta);

WBmwss = gen_bmwss_matrix(N1, N2, alpha, beta, u, v);
WEBmwss = gen_ebmwss_matrix(N1, N2, alpha, beta, u, v);

WGsc = gen_gsc_matrix(N1, N2, alpha, beta, u, v);

W = {{WBmwss}, {WEBmwss},  {WGsc}, {W1DFT, W2DFT}, {W1Golay, W2Golay}};


montNum = 1e5;
frameLen = 1e3;
chan = 'awgn';
if isequal(chan, 'awgn')
    SNR = -12:2:4;
else
    SNR = 0:4:30;
end
ber = zeros(length(W), length(SNR));
for i = 1:length(SNR)
    snr = db2pow(SNR(i));
    for j = 1:length(W)
        w = W{j};
        errorNum = 0;
        parfor k = 1:montNum
            bits = randi([0, 1], 1, frameLen);
            bpsk = 2*bits-1;
            errorNum= errorNum + tranceiver(bpsk, w, snr, alpha, beta, u, v, chan);
        end
        ber(j, i) = errorNum/(montNum*frameLen); 
    end
end

figure
marks = {'-^', '-v', '-o', '-*', '-s'};
for i = 1:size(ber, 1)
    semilogy(SNR, ber(i, :), marks{i}, 'LineWidth', 1.5);
    hold on
end
set(gca, 'LooseInset', [0,0,0,0]);
set(gca, 'FontSize', 12)
legend('BMWSS [12] ', 'EBMWSS [13]', 'Chirp [15]', 'DFT', 'HierBF', fontsize=12)
xlabel('SNR (dB)', FontSize=12)
ylabel('BER', FontSize=12)

