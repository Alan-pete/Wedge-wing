function [beta] = beta_solver(M, gamma, theta)

% solving for shock wave angle (beta)
%delta = 1;          % weak shock

lambda = sqrt((M^2 - 1)^2 - (3 * (1 + ((gamma-1)/2)*M^2) * (1 + ((gamma+1)/2)*M^2) * ((tan((pi/180) * theta))^2)));

chi = ((M^2 - 1)^3 - (9 * (1 + ((gamma-1)/2)*M^2) * (1 + ((gamma-1)/2)*M^2 + ((gamma+1)/4)*M^4) * ((tan((pi/180) * theta))^2))) / lambda^3;

beta_num = (M^2 - 1) + (2 * lambda * cos(((4*pi) + acos(chi)) / 3));
beta_denom = 3 * (1 + ((gamma-1)/2)*M^2) * tan((pi/180) * theta);

beta = (180 / pi) * atan(beta_num / beta_denom);                    % deg

