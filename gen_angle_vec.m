function [theta_vec,phi_vec] = gen_angle_vec(horizon_sample_num, vertical_sample_num, mode, angleLim)
if nargin < 4
    angleLim = [0, 2*pi, -pi/2, pi/2];
end
theta0 = angleLim(1); theta1 = angleLim(2);
phi0 = angleLim(3); phi1 = angleLim(4);
assert(phi1==-phi0, 'Asymmetric elevation!');
if isequal(mode, 'random')
    theta_vec = theta0 + (theta1-theta0)*rand(1, horizon_sample_num);
    phi_vec =phi0 + (phi1-phi0)*rand(1, vertical_sample_num);
%     phi_vec = asin(sin(phi1)*(1-2*rand(1, vertical_sample_num)));
    theta_vec = sort(theta_vec);
    phi_vec = sort(phi_vec);
else
    dTheta = (theta1-theta0)/(horizon_sample_num-1);
    dPhi = (phi1-phi0)/(vertical_sample_num-1);
    theta_vec = theta0:dTheta:theta1;
    phi_vec = phi0:dPhi:phi1;
end
end
