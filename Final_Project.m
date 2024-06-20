clear all; close all; clc;

%% Part 2a: Inviscid L/D Plot and Supporting Analysis
clc;
% given values
delta = 2;                         % half angle (deg)
c = 2;                             % m
Pr = 1;
Rf_lam = sqrt(Pr);
Rf_tur = Pr^(1/3);
gamma = 1.4;

% pulling in the text file for air properties
% Altitude (km)     Pressure (Pa)       T (K)       Rho (kg/m^2)        mu
% (cP)      SonicVel (m/s)
air_matrix = readmatrix('Final Project/US_Standard_Atmosphere.txt');

alpha = linspace(0, 12, 51);
Mach_2a = [1.5, 2, 4, 8];                           % part (a)

% Part (a)
LD_a = zeros(4,length(alpha));

for i = 1:4
    for j = 1:length(alpha)
        % oblique shock across 1-2
        theta_12a = delta + alpha(j);
        [beta_2a, PR_21a, M_2a] = oblique_shock(Mach_2a(i), gamma, theta_12a);
        
        % expansion fan across 2-4
        theta_24a = 2 * delta;
        [M_4a, PR_42a] = expansion_fan(M_2a, gamma, theta_24a);
        PR_41a = PR_21a * PR_42a;
        
        % shock across 1-3
        theta_13a = delta - alpha(j);
        if theta_13a >= 0
            % oblique shock
            [beta_3a, PR_31a, M_3a] = oblique_shock(Mach_2a(i), gamma, abs(theta_13a));
        else
            % expansion fan
            [M_3a, PR_31a] = expansion_fan(Mach_2a(i), gamma, abs(theta_13a));
        end

        % expansion fan across 3-5
        theta_35a = 2 * delta;
        [M_5a, PR_53a] = expansion_fan(M_3a, gamma, theta_35a);
        PR_51a = PR_31a * PR_53a;
        
        % Lift Coefficient
        a = (1 / (2 * cos(pi/180 * delta))) * (1 / ((gamma/2) * Mach_2a(i)^2));
        C_l1a = ((PR_21a - PR_51a) * cos(pi/180 * (delta + alpha(j)))) + ((PR_41a - PR_31a) * cos(pi/180 * (delta - alpha(j))));
        C_la = a * C_l1a;
        
        % Drag Coefficient
        C_d1a = ((PR_21a - PR_51a) * sin(pi/180 * (delta + alpha(j)))) + ((PR_31a - PR_41a) * sin(pi/180 * (delta - alpha(j))));
        C_da = a * C_d1a;
        
        LD_a(i,j) = C_la / C_da;
    end
end

figure(1); 
plot(alpha, LD_a, 'linewidth', 1)
xlabel('Alpha (deg)')
ylabel('Lift/Drag')
title('Inviscid Flow')
legend('M = 1.5', 'M = 2', 'M = 4', 'M = 8')

%% Part 2b: Inviscid L/D Max Plot and Supporting Analysis
clc;
Mach_2b = linspace(1.25, 10, 51);

LD_b = zeros(length(Mach_2b),length(alpha));

% pressure changes with altitude, but for inviscid flow it doesn't affect the lift or drag making them all the same
for i = 1:length(Mach_2b)
    for j = 1:length(alpha)
        % oblique shock across 1-2
        theta_12b = delta + alpha(j);
        [beta_2b, PR_21b, M_2b] = oblique_shock(Mach_2b(i), gamma, theta_12b);
        
        % expansion fan across 2-4
        theta_24b = 2 * delta;
        [M_4b, PR_42b] = expansion_fan(M_2b, gamma, theta_24b);
        PR_41b = PR_21b * PR_42b;
        
        % shock across 1-3
        theta_13b = delta - alpha(j);
        if theta_13b >= 0
            % oblique shock
            [beta_3b, PR_31b, M_3b] = oblique_shock(Mach_2b(i), gamma, abs(theta_13b));
        else
            % expansion fan
            [M_3b, PR_31b] = expansion_fan(Mach_2b(i), gamma, abs(theta_13b));
        end
        
        % expansion fan across 3-5
        theta_35b = 2 * delta;
        [M_5b, PR_53b] = expansion_fan(M_3b, gamma, theta_35b);
        PR_51b = PR_31b * PR_53b;
        
        % Lift Coefficient
        b = (1 / (2 * cos(pi/180 * delta))) * (1 / ((gamma/2) * Mach_2b(i)^2));
        C_l1b = ((PR_21b - PR_51b) * cos(pi/180 * (delta + alpha(j)))) + ((PR_41b - PR_31b) * cos(pi/180 * (delta - alpha(j))));
        C_lb = b * C_l1b;
        
        % Drag Coefficient
        C_d1b = ((PR_21b - PR_51b) * sin(pi/180 * (delta + alpha(j)))) + ((PR_31b - PR_41b) * sin(pi/180 * (delta - alpha(j))));
        C_db = b * C_d1b;
        
        LD_b(i,j) = C_lb / C_db;
    end
end

LD_max10 = zeros(1,length(Mach_2b));
LD_max20 = zeros(1,length(Mach_2b));
LD_max30 = zeros(1,length(Mach_2b));
LD_max40 = zeros(1,length(Mach_2b));
LD_max50 = zeros(1,length(Mach_2b));

for i = 1:length(Mach_2b)
    LD_max10(i) = max(LD_b(i,:));
    LD_max20(i) = max(LD_b(i,:));
    LD_max30(i) = max(LD_b(i,:));
    LD_max40(i) = max(LD_b(i,:));
    LD_max50(i) = max(LD_b(i,:));
end

figure(2);
plot(Mach_2b, LD_max10, Mach_2b, LD_max20, Mach_2b, LD_max30, Mach_2b, LD_max40, Mach_2b, LD_max50, 'linewidth', 1)
xlabel('Mach Number')
ylabel('Lift/Drag Max')
title('Inviscid Flow')
legend('alt = 10 km', 'alt = 20 km', 'alt = 30 km', 'alt = 40 km', 'alt = 50 km')
legend('location', 'southeast')

%% Part 3a - Viscous L/D Plots and Supporting Analysis
clc;
% Text file values for air properties at different altitudes
p_inf = [air_matrix(11,2)/1000, air_matrix(21,2)/1000, air_matrix(31,2)/1000, air_matrix(41,2)/1000, air_matrix(51,2)/1000];
T_inf = [air_matrix(11,3), air_matrix(21,3), air_matrix(31,3), air_matrix(41,3), air_matrix(51,3)];
Mu_inf = [air_matrix(11,5)/1000, air_matrix(21,5)/1000, air_matrix(31,5)/1000, air_matrix(41,5)/1000, air_matrix(51,5)/1000];
V = [air_matrix(11,6), air_matrix(21,6), air_matrix(31,6), air_matrix(41,6), air_matrix(51,6)];
rho_inf = [air_matrix(11,4), air_matrix(21,4), air_matrix(31,4), air_matrix(41,4), air_matrix(51,4)];

Mach_3a = [2, 5, 10];

Re_3 = zeros(5, length(Mach_3a));
T_avg3 = zeros(5, length(Mach_3a));
CD_f_incomp3 = zeros(5, length(Mach_3a));
CD_f_comp3 = zeros(5, length(Mach_3a));

LD_3a = zeros(5, length(Mach_3a), length(alpha));
LD_3a2 = zeros(length(Mach_3a), length(alpha));
LD_3a5 = zeros(length(Mach_3a), length(alpha));
LD_3a10 = zeros(length(Mach_3a), length(alpha));

LD_3a_inv = zeros(length(Mach_3a), length(alpha));
LD_3a2inv = zeros(length(Mach_3a), length(alpha));
LD_3a5inv = zeros(length(Mach_3a), length(alpha));
LD_3a10inv = zeros(length(Mach_3a), length(alpha));

% Part (a)
for i = 1:5
    for j = 1:length(Mach_3a)
        for k = 1:length(alpha)
            % Solving for Reynolds number
            Re_3(i,j) = (rho_inf(i) * (Mach_3a(j) * V(i)) * (c * (cos(pi/180 * alpha(k)) / cos(pi/180 * delta)))) / Mu_inf(i);

            % Solving for average temperature
            T_avg3(i,j) = T_inf(i) * (1 + ((2/9) * (((gamma-1)/2) * Mach_3a(j)^2)));

            % Incompressible skin friction
            CD_f_incomp3(i,j) = 7 / (225 * Re_3(i,j)^(1/7));

            % Compressibility correction (1/7 power rule)
            CD_f_comp3(i,j) = CD_f_incomp3(i,j) / ((T_inf(i)/T_avg3(i,j))^(5/2) * ((T_avg3(i,j) + 120)/(T_inf(i) + 120)))^(1/7);

            % L/D analysis
            % oblique shock across 1-2
            theta_12 = delta + alpha(k);
            [beta_2, PR_21, M_2] = oblique_shock(Mach_3a(j), gamma, theta_12);
            
            % expansion fan across 2-4
            theta_24 = 2 * delta;
            [M_4, PR_42] = expansion_fan(M_2, gamma, theta_24);
            PR_41 = PR_21 * PR_42;
            
            % shock across 1-3
            theta_13 = delta - alpha(k);
            if theta_13 >= 0
                % oblique shock
                [beta_3, PR_31, M_3] = oblique_shock(Mach_3a(j), gamma, abs(theta_13));
            else
                % expansion fan
                [M_3, PR_31] = expansion_fan(Mach_3a(j), gamma, abs(theta_13));
            end
            
            % expansion fan across 3-5
            theta_35 = 2 * delta;
            [M_5, PR_53] = expansion_fan(M_3, gamma, theta_35);
            PR_51 = PR_31 * PR_53;
            
            % Lift Coefficient
            b = (1 / (2 * cos(pi/180 * delta))) * (1 / ((gamma/2) * Mach_3a(j)^2));
            C_l1 = ((PR_21 - PR_51) * cos(pi/180 * (delta + alpha(k)))) + ((PR_41 - PR_31) * cos(pi/180 * (delta - alpha(k))));
            C_l = b * C_l1;
            
            % Drag Coefficient
            C_d1 = ((PR_21 - PR_51) * sin(pi/180 * (delta + alpha(k)))) + ((PR_31 - PR_41) * sin(pi/180 * (delta - alpha(k))));
            C_d = (b * C_d1) + CD_f_comp3(i,j);
            C_dinv = (b * C_d1);
            
            LD_3a(i,j,k) = C_l / C_d;
            LD_3a_inv(i,j,k) = C_l / C_dinv;

            LD_3a2(i,k) = LD_3a(i, 1, k);
            LD_3a5(i,k) = LD_3a(i, 2, k);
            LD_3a10(i,k) = LD_3a(i, 3, k);

            LD_3a2inv(i,k) = LD_3a_inv(i, 1, k);
            LD_3a5inv(i,k) = LD_3a_inv(i, 2, k);
            LD_3a10inv(i,k) = LD_3a_inv(i, 3, k);

        end
    end
end

figure(3);
plot(alpha, LD_3a2, alpha, LD_3a2inv)
xlabel('Alpha (deg)')
ylabel('Lift/Drag')
title('Mach = 2')
legend('alt = 10 km', 'alt = 20 km', 'alt = 30 km', 'alt = 40 km', 'alt = 50 km', 'inviscid')

figure(4);
plot(alpha, LD_3a5, alpha, LD_3a5inv)
xlabel('Alpha (deg)')
ylabel('Lift/Drag')
title('Mach = 5')
legend('alt = 10 km', 'alt = 20 km', 'alt = 30 km', 'alt = 40 km', 'alt = 50 km', 'inviscid')

figure(5);
plot(alpha, LD_3a10, alpha, LD_3a10inv)
xlabel('Alpha (deg)')
ylabel('Lift/Drag')
title('Mach = 10')
legend('alt = 10 km', 'alt = 20 km', 'alt = 30 km', 'alt = 40 km', 'alt = 50 km', 'inviscid')

%% Part 3b - Viscous L/D Max Plots and Supporting Analysis
clc;

Mach_3b = linspace(1.25, 10, 51);

LD_3b = zeros(5, length(Mach_3b), length(alpha));
LD_3b10 = zeros(length(Mach_3b), length(alpha));
LD_3b20 = zeros(length(Mach_3b), length(alpha));
LD_3b30 = zeros(length(Mach_3b), length(alpha));
LD_3b40 = zeros(length(Mach_3b), length(alpha));
LD_3b50 = zeros(length(Mach_3b), length(alpha));

for i = 1:5
    for j = 1:length(Mach_3b)
        for k = 1:length(alpha)
            % Solving for Reynolds number
            Re_3(i,j) = (rho_inf(i) * (Mach_3b(j) * V(i)) * (c * (cos(pi/180 * alpha(k)) / cos(pi/180 * delta)))) / Mu_inf(i);

            % Solving for average temperature
            T_avg3(i,j) = T_inf(i) * (1 + ((2/9) * (((gamma-1)/2) * Mach_3b(j)^2)));

            % Incompressible skin friction
            CD_f_incomp3(i,j) = 7 / (225 * Re_3(i,j)^(1/7));

            % Compressibility correction (1/7 power rule)
            CD_f_comp3(i,j) = CD_f_incomp3(i,j) / ((T_inf(i)/T_avg3(i,j))^(5/2) * ((T_avg3(i,j) + 120)/(T_inf(i) + 120)))^(1/7);

            % L/D analysis
            % oblique shock across 1-2
            theta_12 = delta + alpha(k);
            [beta_2, PR_21, M_2] = oblique_shock(Mach_3b(j), gamma, theta_12);
            
            % expansion fan across 2-4
            theta_24 = 2 * delta;
            [M_4, PR_42] = expansion_fan(M_2, gamma, theta_24);
            PR_41 = PR_21 * PR_42;
            
            % shock across 1-3
            theta_13 = delta - alpha(k);
            if theta_13 >= 0
                % oblique shock
                [beta_3, PR_31, M_3] = oblique_shock(Mach_3b(j), gamma, abs(theta_13));
            else
                % expansion fan
                [M_3, PR_31] = expansion_fan(Mach_3b(j), gamma, abs(theta_13));
            end
            
            % expansion fan across 3-5
            theta_35 = 2 * delta;
            [M_5, PR_53] = expansion_fan(M_3, gamma, theta_35);
            PR_51 = PR_31 * PR_53;
            
            % Lift Coefficient
            b = (1 / (2 * cos(pi/180 * delta))) * (1 / ((gamma/2) * Mach_3b(j)^2));
            C_l1 = ((PR_21 - PR_51) * cos(pi/180 * (delta + alpha(k)))) + ((PR_41 - PR_31) * cos(pi/180 * (delta - alpha(k))));
            C_l = b * C_l1;
            
            % Drag Coefficient
            C_d1 = ((PR_21 - PR_51) * sin(pi/180 * (delta + alpha(k)))) + ((PR_31 - PR_41) * sin(pi/180 * (delta - alpha(k))));
            C_d = (b * C_d1) + CD_f_comp3(i,j);
            C_dinv = (b * C_d1);
            
            LD_3b(i,j,k) = C_l / C_d;

            LD_3b10(j,k) = LD_3b(1, j, k);
            LD_3b20(j,k) = LD_3b(2, j, k);
            LD_3b30(j,k) = LD_3b(3, j, k);
            LD_3b40(j,k) = LD_3b(4, j, k);
            LD_3b50(j,k) = LD_3b(5, j, k);

        end
    end
end

LD_max10 = zeros(1,length(Mach_3b));
LD_max20 = zeros(1,length(Mach_3b));
LD_max30 = zeros(1,length(Mach_3b));
LD_max40 = zeros(1,length(Mach_3b));
LD_max50 = zeros(1,length(Mach_3b));

for i = 1:length(Mach_2b)
    LD_max10(i) = max(LD_3b10(i,:));
    LD_max20(i) = max(LD_3b20(i,:));
    LD_max30(i) = max(LD_3b30(i,:));
    LD_max40(i) = max(LD_3b40(i,:));
    LD_max50(i) = max(LD_3b50(i,:));
end

figure(6);
plot(Mach_3b, LD_max10, Mach_3b, LD_max20, Mach_3b, LD_max30, Mach_3b, LD_max40, Mach_3b, LD_max50)
xlabel('Mach Number')
ylabel('Lift/Drag Max')
title('Viscous L/D Max')
legend('alt = 10 km', 'alt = 20 km', 'alt = 30 km', 'alt = 40 km', 'alt = 50 km')

%% Part 4a: Transitional Flow L/D Plot and Supporting Analysis

clc;

Mach_4a = linspace(1.25, 10, 51);

Re_4 = zeros(5, length(Mach_4a));
T_avg4 = zeros(5, length(Mach_4a));
CD_f_incomp4 = zeros(5, length(Mach_4a));
CD_f_comp4 = zeros(5, length(Mach_4a));

LD_4b = zeros(5, length(Mach_4a), length(alpha));
LD_4b10 = zeros(length(Mach_4a), length(alpha));
LD_4b20 = zeros(length(Mach_4a), length(alpha));
LD_4b30 = zeros(length(Mach_4a), length(alpha));
LD_4b40 = zeros(length(Mach_4a), length(alpha));
LD_4b50 = zeros(length(Mach_4a), length(alpha));

for i = 1:5
    for j = 1:length(Mach_4a)
        for k = 1:length(alpha)
            % Solving for Reynolds number
            Re_4(i,j) = (rho_inf(i) * (Mach_4a(j) * V(i)) * (c * (cos(pi/180 * alpha(k)) / cos(pi/180 * delta)))) / Mu_inf(i);

            % Laminar compressibility correction
            if Re_4(i,j) < 500000
                % Solving for average temperature
                T_avg4(i,j) = T_inf(i) * (1 + ((16/5) * (((gamma-1)/2) * Mach_4a(j)^2)));

                % Incompressible skin friction
                CD_f_incomp4(i,j) = 1.328 / (Re_4(i,j)^(1/2));

                % Compressibility correction
                CD_f_comp4(i,j) = CD_f_incomp4(i,j) / ((T_inf(i)/T_avg4(i,j))^(5/2) * ((T_avg4(i,j) + 120)/(T_inf(i) + 120)))^(1/2);
            end

            % Turbulent and transitional compressibility correction
            if Re_4(i,j) >= 500000
                % Solving for average temperature
                T_avg4(i,j) = T_inf(i) * (1 + ((2/9) * (((gamma-1)/2) * Mach_4a(j)^2)));

                % Incompressible skin friction
                CD_f_incomp4(i,j) = 7 / (225 * Re_4(i,j)^(1/7));
    
                % Compressibility correction (1/7 power rule)
                CD_f_comp4(i,j) = CD_f_incomp4(i,j) / ((T_inf(i)/T_avg4(i,j))^(5/2) * ((T_avg4(i,j) + 120)/(T_inf(i) + 120)))^(1/7);
            end

            % L/D analysis
            % oblique shock across 1-2
            theta_12 = delta + alpha(k);
            [beta_2, PR_21, M_2] = oblique_shock(Mach_4a(j), gamma, theta_12);
            
            % expansion fan across 2-4
            theta_24 = 2 * delta;
            [M_4, PR_42] = expansion_fan(M_2, gamma, theta_24);
            PR_41 = PR_21 * PR_42;
            
            % shock across 1-3
            theta_13 = delta - alpha(k);
            if theta_13 >= 0
                % oblique shock
                [beta_3, PR_31, M_3] = oblique_shock(Mach_4a(j), gamma, abs(theta_13));
            else
                % expansion fan
                [M_3, PR_31] = expansion_fan(Mach_4a(j), gamma, abs(theta_13));
            end
            
            % expansion fan across 3-5
            theta_35 = 2 * delta;
            [M_5, PR_53] = expansion_fan(M_3, gamma, theta_35);
            PR_51 = PR_31 * PR_53;
            
            % Lift Coefficient
            b = (1 / (2 * cos(pi/180 * delta))) * (1 / ((gamma/2) * Mach_4a(j)^2));
            C_l1 = ((PR_21 - PR_51) * cos(pi/180 * (delta + alpha(k)))) + ((PR_41 - PR_31) * cos(pi/180 * (delta - alpha(k))));
            C_l = b * C_l1;
            
            % Drag Coefficient
            C_d1 = ((PR_21 - PR_51) * sin(pi/180 * (delta + alpha(k)))) + ((PR_31 - PR_41) * sin(pi/180 * (delta - alpha(k))));
            C_d = (b * C_d1) + CD_f_comp3(i,j);
            C_dinv = (b * C_d1);
            
            LD_4b(i,j,k) = C_l / C_d;

            LD_4a10(j,k) = LD_4b(1, j, k);
            LD_4a20(j,k) = LD_4b(2, j, k);
            LD_4a30(j,k) = LD_4b(3, j, k);
            LD_4a40(j,k) = LD_4b(4, j, k);
            LD_4a50(j,k) = LD_4b(5, j, k);

        end
    end
end

LD_max10 = zeros(1,length(Mach_4a));
LD_max20 = zeros(1,length(Mach_4a));
LD_max30 = zeros(1,length(Mach_4a));
LD_max40 = zeros(1,length(Mach_4a));
LD_max50 = zeros(1,length(Mach_4a));

for i = 1:length(Mach_4a)
    LD_max10(i) = max(LD_4a10(i,:));
    LD_max20(i) = max(LD_4a20(i,:));
    LD_max30(i) = max(LD_4a30(i,:));
    LD_max40(i) = max(LD_4a40(i,:));
    LD_max50(i) = max(LD_4a50(i,:));
end

figure(7);
plot(Mach_4a, LD_max10, Mach_4a, LD_max20, Mach_4a, LD_max30, Mach_4a, LD_max40, Mach_4a, LD_max50)
xlabel('Mach Number')
ylabel('Lift/Drag Max')
title('Viscous L/D Max')
legend('alt = 10 km', 'alt = 20 km', 'alt = 30 km', 'alt = 40 km', 'alt = 50 km')

% Part (b)

figure(8);
plot(Re_4(1,:), LD_max10, Re_4(2,:), LD_max20, Re_4(3,:), LD_max30, Re_4(4,:), LD_max40, Re_4(5,:), LD_max50)
xline(500000, '-','Laminar-Transition')
xline(10^7, '-', 'Transition-Turbulent')
set(gca, 'XScale', 'log')
xlabel('Reynolds Number Re')
ylabel('Lift/Drag Max')
title('Viscous L/D Max')
legend('alt = 10 km', 'alt = 20 km', 'alt = 30 km', 'alt = 40 km', 'alt = 50 km')






























