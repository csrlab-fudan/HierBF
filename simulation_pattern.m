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

WBmwss = gen_bmwss_matrix(N1, N2, alpha, beta, u, v);
WEBmwss = gen_ebmwss_matrix(N1, N2, alpha, beta, u, v);

WGsc = gen_gsc_matrix(N1, N2, alpha, beta, u, v);

W = {{WBmwss}, {WEBmwss},  {WGsc}, {W1Golay, W2Golay}};

%% 
u = 0; v = 0;
angleLim = [pi/3, 2*pi/3, -pi/6, pi/6];
beamSampleHorizonNum = 180;
beamSampleVerticalNum = 180;
beamSampleLength = beamSampleHorizonNum*beamSampleVerticalNum;
broadbeampattern = zeros(beamSampleHorizonNum, beamSampleVerticalNum);
mode = 'uniform'; % the mode of choosing discrete angles
[beamThetaVec, beamPhiVec] = gen_angle_vec(beamSampleHorizonNum, beamSampleVerticalNum, mode, angleLim);
[~, indexTargetTheta] = min(abs(targetTheta-beamThetaVec));
[~, indexTargetPhi] = min(abs(targetPhi-beamPhiVec));
maskTheta = zeros(size(beamThetaVec));
maskPhi = zeros(size(beamPhiVec));
% stepTheta = round(1/2*abs(acos(u-alpha)-acos(u+alpha))/(2*pi)*length(beamThetaVec));
% stepPhi = round(1/2*abs(asin(v-beta)-asin(v+beta))/pi*length(beamPhiVec));
idxTheta = beamThetaVec <= acos(u-alpha) & beamThetaVec >=acos(u+alpha);
idxPhi = beamPhiVec <= asin(v+beta) & beamPhiVec >= asin(v-beta);
maskTheta(idxTheta) = 1/alpha/beta;
maskPhi(idxPhi) = 1/alpha/beta;
[~, right] = max(beamThetaVec >= acos(u-alpha));
[~, left] = max(beamThetaVec >= acos(u+alpha));
[~, top] = max(beamPhiVec >= asin(v-beta));
[~, down] = max(beamPhiVec >= asin(v+beta));
boundx = [left, left, right, right, left];
boundy = [down, top, top, down, down];

close all
titles = {'BMWSS [13] ', 'EBMWSS [14]', 'Chirp [16]', 'HierBF', 'Mask'};
f1 = figure;
t1 = tiledlayout('flow');
f2 = figure;
t2 = tiledlayout('flow');
f3 = figure;
t3 = tiledlayout('flow');
f4 = figure;
t4 = tiledlayout('flow');
for ii = 1:length(W)
    w = W{ii};
    if length(w) == 2
        for i = 1:beamSampleHorizonNum
            temp = cos(beamThetaVec(i));
            w1 = w{1}; w2 = w{2};
            for j = 1:beamSampleVerticalNum
                u = sin(beamPhiVec(j));
                v = temp*cos(beamPhiVec(j));
                F  = exp(-1j*pi*(u*(0:N1-1).' + (v*(0:N2-1))));
                W1 = w1; W2 = w2;
                broadbeampattern(i, j) = 0.5*(abs(F(:).'*W1(:))^2 + abs(F(:).'*W2(:))^2);
            end
        end
    else
        for i = 1:beamSampleHorizonNum
            temp = cos(beamThetaVec(i));
            w1 = w{1};
            for j = 1:beamSampleVerticalNum
                u = sin(beamPhiVec(j));
                v = temp*cos(beamPhiVec(j));
                F  = exp(-1j*pi*(u*(0:N1-1).' + (v*(0:N2-1))));
                W1 = w1;
                broadbeampattern(i, j) = abs(F(:).'*W1(:))^2;
            end
        end
    end
%     figure
%     feature('DefaultCharacterSet','UTF-8');
%     [phi, theta] = meshgrid(beamPhiVec, beamThetaVec);
%     [x,y,z] = sph2cart(theta,phi,broadbeampattern);
%     p = mesh(x, y, z, broadbeampattern);
% %     view(180, 90)
%     axis equal
%     xlim([0, inf])
%     ylim([0, inf])
%     ax = gca;
%     set(gca,'XTick',[])
%     set(gca,'YTick',[])
%     set(gca,'ZTick',[])
%     set(gca, 'LooseInset', [0,0,0,0]);
%     set(gcf, "Position", [0, 0, 600, 600]);
%     hidden off;
%     light('position',[0 0 -5]);
%     hold on

    %  cross-section in azimuth plane

    figure(f1)
%     [x, y] = pol2cart(beamThetaVec, broadbeampattern(:, indexTargetPhi).');
%     plot(x, y)
%     axis equal
%     box off; axis off
%     set(gca, 'LooseInset', [0,0,0,0]);
%     set(gcf, "Position", [0, 0, 900, 300]);
%     view(90, 90)
    nexttile
    plot(180/pi*beamThetaVec, broadbeampattern(:, indexTargetPhi).', LineWidth=1.0);
    hold on
    plot(180/pi*beamThetaVec, maskTheta, '--', LineWidth=1.0);
    title(titles{ii})




    % cross-section in elevation plane
    figure(f2)
%     [x, y] = pol2cart(beamPhiVec, broadbeampattern(indexTargetTheta, :));
%     plot(x, y)
%     axis equal
%     box off; axis off
%     set(gca, 'LooseInset', [0,0,0,0]);
%     set(gcf, "Position", [0, 0, 900, 300]);
    nexttile
    plot(180/pi*beamPhiVec, broadbeampattern(indexTargetTheta, :), LineWidth=1.0);
    hold on
    plot(180/pi*beamPhiVec, maskPhi, '--', LineWidth=1.0);
    title(titles(ii))

    figure(f3)
    nexttile
    imagesc(broadbeampattern)
    xlabel('$\varphi / \circ$', 'Interpreter','latex')
    ylabel('$\theta / \circ$', 'Interpreter','latex')
    xtickslabels = linspace(angleLim(3)*180/pi, angleLim(4)*180/pi, 4);
    xticks = linspace(1, size(broadbeampattern, 2), numel(xtickslabels));
    set(gca, 'XTick', xticks, 'XTickLabel', xtickslabels)
    ytickslabels = linspace(angleLim(1)*180/pi, angleLim(2)*180/pi, 4);
    yticks = linspace(1, size(broadbeampattern, 1), numel(ytickslabels));
    set(gca, 'YTick', yticks, 'YTickLabel', ytickslabels)
    hold on
    plot(boundx, boundy, '--r', LineWidth=1.0)
    colorbar
    title(titles(ii))

    figure(f4)
    nexttile
    plot_phase(W1)
    title(titles(ii), FontSize=9)
end
figure(f1)
t1.TileSpacing = 'compact';
t1.Padding = 'compact';
xlabel(t1, '$\theta / \circ$', 'Interpreter','latex');
ylabel(t1, 'Power')
figure(f2)
t2.TileSpacing = 'compact';
t2.Padding = 'compact';
xlabel(t2, '$\varphi$', 'Interpreter','latex');
ylabel(t2, 'Power')
% figure(f4)
% t4.TileSpacing = 'compact';
% t4.Padding = 'compact';

function plot_phase(A)
a = angle(A(:));
polarscatter(a, ones(length(a), 1), 20);
rlim([0, 1])
ax = gca;
ax.RTickLabel = {};
ax.RGrid = 'off';
% ax.ThetaGrid = 'off';
set(gca, "FontSize", 9)
% set(gca,'LooseInset',get(gca,'TightInset'))
end

