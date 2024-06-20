function [M] = Prandtl_Meyer(M_1, gamma, alpha, n)

% determining the first guess for Mach number
m_init = M_1 + 0.1;

% starting a list for the iterations of the Mach number, starts with the
% initial guess
M_n = zeros(1, n);
M_n(1) = m_init;

% V(M_1)
vm1_1 = sqrt((gamma + 1) / (gamma - 1));
vm1_2 = atan(sqrt(((gamma - 1) / (gamma + 1)) * (M_1^2 - 1)));
vm1_3 = atan(sqrt(M_1^2 - 1));

vm_1 = (vm1_1 * vm1_2) - vm1_3;

% V(M_2)
vm_2 = (alpha * (pi / 180)) + vm_1;

error = 1;

% for loop to iterate through
while error > 0.001
    for i = 2:n
        % V(M_2) j
        vm2_1j = sqrt((gamma + 1) / (gamma - 1));
        vm2_2j = atan(sqrt(((gamma - 1) / (gamma + 1)) * (M_n(i-1)^2 - 1)));
        vm2_3j = atan(sqrt(M_n(i-1)^2 - 1));

        vm_2j = (vm2_1j * vm2_2j) - vm2_3j;

        % Prandtl-Meyer derivative
        dv_dm = (1 / M_n(i-1)) * ((sqrt(M_n(i-1)^2 - 1)) / (1 + (((gamma - 1) / 2) * M_n(i-1)^2)));

        % full Newton solver
        M_n(i) = ((vm_2 - vm_2j) / dv_dm) + M_n(i-1);

        % error to decide when to drop from loop
        error = abs((M_n(i) - M_n(i-1)) / M_n(i));

    end
end

M = M_n(n);