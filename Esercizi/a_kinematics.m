close all; clearvars; clc;

% Load the URDF model
robot = importrobot('RPR_zyx_LucaPonti.urdf');

% home configuration
config = homeConfiguration(robot);
% figure
% show(robot, config);

config(1).JointPosition = 0;
config(2).JointPosition = 0;
config(3).JointPosition = pi/2;
q1 = config(1).JointPosition;
q2 = config(2).JointPosition;
q3 = config(3).JointPosition;

% figure
show(robot, config);


disp("====== Direct Kinematic ======");

lb = 0.15;
l1 = 0.4;
l3 = 0.24;

% [theta alfa a d]
syms theta1 theta2 theta3 real
DH_matrix = [0 0 0 lb;
             theta1 pi / 2 l1 0;
             pi / 2 -pi / 2 0 theta2 + 0.3;
             theta3 - pi / 2 0 l3 0;
             pi / 2 pi / 2 0 0];

syms a alpha d theta real
A = [cos(theta), -sin(theta) * cos(alpha), sin(theta) * sin(alpha), a * cos(theta);
     sin(theta), cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta);
     0, sin(alpha), cos(alpha), d;
     0, 0, 0, 1];

A_b_0 = subs(A, [theta alpha a d], DH_matrix(1, :));
A_0_1 = subs(A, [theta alpha a d], DH_matrix(2, :));
A_b_1 = A_b_0 * A_0_1;
A_1_2 = subs(A, [theta alpha a d], DH_matrix(3, :));
A_b_2 = A_b_0 * A_0_1 * A_1_2;
A_2_3 = subs(A, [theta alpha a d], DH_matrix(4, :));
A_b_3 = A_b_0 * A_0_1 * A_1_2 * A_2_3;
A_3_ee = subs(A, [theta alpha a d], DH_matrix(5, :));
A_b_ee = A_b_0 * A_0_1 * A_1_2 * A_2_3 * A_3_ee;

% just for the visualization of the DH table
Frames = {'Base --> 0'; "0 --> 1"; "1 --> 2"; "2 --> 3"; '3 --> End-effector'};
DH_cell_matrix = [Frames, num2cell(DH_matrix)];
DH_table = cell2table(DH_cell_matrix, 'VariableNames', {'Frames', 'Theta', 'Alpha', 'a', 'd'});
disp("DH table:")
disp(DH_table);

disp("====== Inverse Kinematic ======");
% when you see the underscore after a variable, it refers to the internal variables of the inverse kinematics

% this works with the position of the ee equal to [0.4; -0.54; 0.15]. The expected values for the joint variables are [0; 0; 0].
th2_ = 0; % found by hand, see exercise on paper
th1_ = solve(0.4 * cos(theta1) + 0.3 * sin(theta1) + th2_ * sin(theta1) == 0.4, theta1, theta1 < pi / 2); % it's 0, as exprected
A_b_ee_ = subs(A_b_ee, [theta1 theta2], [th1_ th2_]); % it's 0, as exprected
th3_ = solve(A_b_ee_(2, 4) == -0.54); % it's 0, as exprected

fprintf('Joint variables founded with Inverse kinematics: [%.2f, %.2f, %.2f]\n', th1_, th2_, th3_);

disp("====== Geometric Jacobian ======");

z1 = A_b_0(1:3, 1:3) * [0; 0; 1];
z2 = A_b_1(1:3, 1:3) * [0; 0; 1];
z3 = A_b_2(1:3, 1:3) * [0; 0; 1];

Pee_P1 = A_b_ee(1:3, 4) - A_b_0(1:3, 4);
Pee_P3 = A_b_ee(1:3, 4) - A_b_2(1:3, 4);

J_geometric = simplify([cross(z1, Pee_P1), z2, cross(z3, Pee_P3);
                        z1, zeros(3, 1), z3]); % symbolic geometric jacobian

disp(double(subs(J_geometric, [theta1, theta2, theta3], [q1 q2 q3])));
J_geometric_toolbox = geometricJacobian(robot, config, "ee");

disp("Found by hand:"); % it's equal to the one found with the toolbox!
disp(J_geometric);
disp("Found with the toolbox:");
disp(J_geometric_toolbox);

disp("====== Analytical Jacobian ======");

Px_ee = A_b_ee(1, 4);
Py_ee = A_b_ee(2, 4);
Pz_ee = A_b_ee(3, 4);

% J_analitycal = [diff(Px_ee, theta1), diff(Px_ee, theta2), diff(Px_ee, theta3);
%                 diff(Py_ee, theta1), diff(Py_ee, theta2), diff(Py_ee, theta3);
%                 diff(Pz_ee, theta1), diff(Pz_ee, theta2), diff(Pz_ee, theta3);
%                 ];

J_pos = simplify(jacobian([Px_ee; Py_ee; Pz_ee], [theta1, theta2, theta3])); % equal to the corresponding part in the geometric jacobian

phi = atan2(A_b_ee(2, 3), A_b_ee(1, 3));
gamma = acos(A_b_ee(3, 3));
psi = atan2(A_b_ee(3, 1), A_b_ee(3, 2));
J_ang = simplify(jacobian([phi, gamma, psi], [theta1, theta2, theta3]));

J_analytical_tot = [J_pos; J_ang]; % symbolic analytical jacobian
J_analytical_tot = subs(J_analytical_tot, [theta1 theta2 theta3], [q1 q2 q3]);
disp(double(J_analytical_tot));
