function [M_2, PR] = expansion_fan(M_1, gamma, theta)

M_2 = Prandtl_Meyer(M_1, gamma, theta, 20);

PR_up = 1 + (((gamma - 1) / 2) * M_1^2);
PR_down = 1 + (((gamma - 1) / 2) * M_2^2);

PR = (PR_up / PR_down)^(gamma / (gamma - 1));