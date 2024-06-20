function [beta, PR, M_2] = oblique_shock(M_1, gamma, theta)

% shock angle
beta = beta_solver(M_1, gamma, theta);

% Normal Mach component 1
M_n1 = M_1 * sin(beta * pi/180);

% Temperature ratio 1-2
TR_1 = 1 + (((2*gamma)/(gamma+1)) * ((M_1 * sin(beta * pi/180))^2 - 1));
TR_2 = (2 + ((gamma-1) * (M_1 * sin(beta * pi/180))^2)) / ((gamma+1) * (M_1 * sin(beta * pi/180))^2);

TR = TR_1 * TR_2;

% Normal Mach 2 component
M_n2_num = 1 + ((gamma-1)/2)*(M_n1^2);
M_n2_den = (gamma * (M_n1^2)) - ((gamma-1)/2);

M_n2 = sqrt(M_n2_num / M_n2_den);

% Tangential Mach 2 component
M_t2 = M_1 * cos(beta * pi/180) * sqrt(1/TR);

% Total Mach behind shock
M_2 = sqrt(M_n2^2 + M_t2^2);

% pressure ratio behind shock
PR = 1 + (((2*gamma)/(gamma+1)) * (M_n1^2 - 1));
