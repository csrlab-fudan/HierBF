
clear;

antRow = 40;
antCol = 36;
targetTheta = pi/4;
targetPhi = 0;
u0 = cos(targetPhi)*cos(targetTheta);
v0 = sin(targetPhi);
Rs = [40, 10, 4, 1];
Cs = [36, 12, 4, 1];
close all
f = figure;
t = tiledlayout('flow');
titles = ["(a)", "(b)", "(c)", "(d)"];
for ii = 1:length(Rs)
    R = Rs(ii); C = Cs(ii); % virtual size
    beamSampleHorizonNum = 900;
    beamSampleVerticalNum = 450;
    beamSampleLength = beamSampleHorizonNum*beamSampleVerticalNum;
    broadbeampattern = zeros(beamSampleHorizonNum, beamSampleVerticalNum);
    mode = 'uniform'; % the mode of choosing discrete angles
    [beamThetaVec, beamPhiVec] = gen_angle_vec(beamSampleHorizonNum, beamSampleVerticalNum, mode);
    
    [W1, W2] = gen_golay_hier(antRow, antCol, R, C, v0, u0);
    
    for i = 1:beamSampleHorizonNum
        temp = cos(beamThetaVec(i));
        parfor j = 1:beamSampleVerticalNum
            v = sin(beamPhiVec(j));
            u = temp*cos(beamPhiVec(j));
            F  = exp(-1j*pi*(v*(0:antRow-1).' + (u*(0:antCol-1))));
            broadbeampattern(i, j) = 1/2*(abs(F(:).'*W1(:))^2 + abs(F(:).'*W2(:))^2);
        end
    end
    

    feature('DefaultCharacterSet','UTF-8');
    [phi, theta] = meshgrid(beamPhiVec, beamThetaVec);
    [x,y,z] = sph2cart(theta,phi,broadbeampattern);
    nexttile
    p = mesh(x, y, z, broadbeampattern);
    axis off
    if ii==4
        caxis([0.9, 1.1])
    end
    colorbar
%     xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');
    view(180, 90)
    axis equal
    xlim([0, inf])
    ylim([0, inf])
    ax = gca;
%     set(ax, 'FontSize', 16)
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'ZTick',[])
%     set(gcf, "Position", [0, 0, 380, 300]);
    
%     ax.XAxis.Exponent = ii+1;
%     ax.YAxis.Exponent = ii+1;
%     hidden off;
    light('position',[0 0 -5]);
    text(max(x(:))/2, max(y(:))*1.1, titles(ii));
    hold on
end
t.TileSpacing = 'compact';
t.Padding = 'compact';
