function truss_analysis()

    % Given data
    AB = 350;
    BC = 346;
    AC = 554;
    theta_C = deg2rad(37.24);
    theta_B = deg2rad(106);
    theta_A = deg2rad(36.76);
    EA_AB = 8;
    EA_AC = 10;
    EA_BC = 10;

    % Joint coordinates
    A = [0, 0];
    B = [AB, 0];
    C = [AB + AC * cos(theta_C), AC * sin(theta_C)];

    % Element lengths and angles
    L_AB = AB;
    L_BC = BC;
    L_AC = AC;

    % Stiffness matrices for each element
    K_AB = truss_stiffness(EA_AB, L_AB, theta_A);
    K_BC = truss_stiffness(EA_BC, L_BC, theta_B);
    K_AC = truss_stiffness(EA_AC, L_AC, theta_C);

    % Global stiffness matrix assembly
    K_global = zeros(4);

    % Add contributions from each element
    K_global(1:2, 1:2) = K_global(1:2, 1:2) + K_AB(1:2, 1:2);
    K_global(3:4, 3:4) = K_global(3:4, 3:4) + K_BC(1:2, 1:2);
    K_global(1:2, 3:4) = K_global(1:2, 3:4) - K_AB(1:2, 3:4);
    K_global(3:4, 1:2) = K_global(3:4, 1:2) - K_AB(3:4, 1:2);
    K_global(3:4, 3:4) = K_global(3:4, 3:4) + K_AC(1:2, 1:2);

    % Apply boundary conditions
    K_reduced = K_global(3:4, 3:4);
    F = [0; 0; 1; 1]; % External forces at joint A
    U_reduced = K_reduced \ F(3:4);

    % Assemble full displacement vector
    U_full = [0; 0; U_reduced(1); U_reduced(2)];

    % Calculate reactions
    F_reactions = K_global(1:2, :) * U_full;

    % Display results
    disp('Displacements (U):');
    disp(U_full);
    disp('Reactions (F):');
    disp(F_reactions);

end

function K = truss_stiffness(EA, L, theta)
    c = cos(theta);
    s = sin(theta);

    % Stiffness matrix for a 2D truss element
    K = (EA / L) * [
        c^2, c * s, -c^2, -c * s;
        c * s, s^2, -c * s, -s^2;
        -c^2, -c * s, c^2, c * s;
        -c * s, -s^2, c * s, s^2
    ];

    disp('Stiffness Matrix:');
    disp(K);
end
