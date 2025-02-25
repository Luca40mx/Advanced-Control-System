clearvars; close all; clc;

addpath("utils\");

%% initial values and constants
jointValues = [pi / 2 -0.2 pi / 3];
velocityValues = [5 2 8];
accelerationValues = [4 1 2];

% link lengths
l_base = 0.15;
l1 = 0.4;
l2 = 0.3;
l3 = 0.24;

%% robot model
robot = importrobot("RPR_zyx_LucaPonti.urdf", "urdf", DataFormat = "row");

% home configuration
config = robot.homeConfiguration;
config(1) = jointValues(1);
config(2) = jointValues(2);
config(3) = jointValues(3);

%% forward kinematics

% [theta alfa a d]
syms theta1 theta2 theta3 real
syms d_theta1 d_theta2 d_theta3 real
syms dd_theta1 dd_theta2 dd_theta3 real
syms th1(t) th2(t) th3(t)
syms d_th1(t) d_th2(t) d_th3(t)
syms dd_th1(t) dd_th2(t) dd_th3(t)

DH_matrix = [0 0 0 l_base;
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

%% Geometric Jacobian

z0 = A_b_0(1:3, 3);
z1 = A_b_1(1:3, 3);
z2 = A_b_2(1:3, 3);

Pee_P1 = A_b_ee(1:3, 4) - A_b_0(1:3, 4);
Pee_P3 = A_b_ee(1:3, 4) - A_b_2(1:3, 4);

J_geometric = simplify([cross(z0, Pee_P1), z1, cross(z2, Pee_P3);
                        z0, zeros(3, 1), z2]); % symbolic geometric jacobian

%% Analytical Jacobian
Px_ee = A_b_ee(1, 4);
Py_ee = A_b_ee(2, 4);
Pz_ee = A_b_ee(3, 4);

J_pos = simplify(jacobian([Px_ee; Py_ee; Pz_ee], [theta1, theta2, theta3])); % equal to the corresponding part in the geometric jacobian

phi = atan2(A_b_ee(2, 3), A_b_ee(1, 3));
gamma = acos(A_b_ee(3, 3));
psi = atan2(A_b_ee(3, 2), -A_b_ee(3, 1));
J_ang = simplify(jacobian([phi, gamma, psi], [theta1, theta2, theta3]));

J_analytical = simplify([J_pos; J_ang]); % symbolic analytical jacobian

J_analytical_transpose = J_analytical.';

J_analytical_inverse = pinv(J_analytical);

J_analytical_time_dependent = subs(J_analytical, [theta1 theta2 theta3], [th1(t) th2(t) th3(t)]);
J_analytical_dot_time_dependent = diff(J_analytical_time_dependent, t);

J_analitycal_dot = subs(J_analytical_dot_time_dependent, [diff(th1(t), t), diff(th2(t), t), diff(th3(t), t)], [d_theta1, d_theta2, d_theta3]); %#ok<*NASGU>
J_analitycal_dot = subs(J_analitycal_dot, [th1(t) th2(t) th3(t)], [theta1 theta2 theta3]);

%% inertia matrix

link1 = robot.getBody("Link1");
mass1 = link1.Mass;
radius1 = 0.02;
[cylinderInertia1, cylinderInertiaMatrix1, traslatedCylinderInertia1] = inertia_cylinder(radius1, 0, l1, mass1);

link2 = robot.getBody("Link2");
mass2 = link2.Mass;
[prismInertia2, prismInertiaMatrix2, traslatedPrismInertia2] = inertia_prism(l2, 0.03, 0.03, mass2);

link3 = robot.getBody("Link3");
mass3 = link3.Mass;
radius3 = 0.02;
[cylinderInertia3, cylinderInertiaMatrix3, traslatedCylinderInertia3] = inertia_cylinder(radius3, 0, l3, mass3);

%% center of mass

com1 = A_b_1 * [-l1 / 2; 0; 0; 1];
com1 = com1(1:3, 1);
com2 = A_b_2 * [0; 0.15; 0; 1];
com2 = com2(1:3, 1);
com3 = A_b_3 * [-l3 / 2; 0; 0; 1];
com3 = com3(1:3, 1);

%% partial jacobian w.r.t. COMs

P0 = A_b_0(1:3, 4);
P1 = A_b_1(1:3, 4);
P2 = A_b_2(1:3, 4);

j_linear_link1 = [cross(z0, com1 - P0), [0; 0; 0], [0; 0; 0]];
j_angular_link1 = [z0, [0; 0; 0], [0; 0; 0]];
j_link1 = [j_linear_link1; j_angular_link1];

j_linear_link2 = [cross(z0, com2 - P0), z1, [0; 0; 0]];
j_angular_link2 = [z0, [0; 0; 0], [0; 0; 0]];
j_link2 = [j_linear_link2; j_angular_link2];

j_linear_link3 = [cross(z0, com3 - P0), z1, cross(z2, com3 - P2)];
j_angular_link3 = [z0, [0; 0; 0], z2];
j_link3 = [j_linear_link3; j_angular_link3];

%% kinetic energy

element1 = mass1 * (j_linear_link1)' * j_linear_link1 + (j_angular_link1)' * A_b_1(1:3, 1:3) * traslatedCylinderInertia1 * (A_b_1(1:3, 1:3))' * j_angular_link1;
element2 = mass2 * (j_linear_link2)' * j_linear_link2 + (j_angular_link2)' * A_b_2(1:3, 1:3) * traslatedPrismInertia2 * (A_b_2(1:3, 1:3))' * j_angular_link2;
element3 = mass3 * (j_linear_link3)' * j_linear_link3 + (j_angular_link3)' * A_b_3(1:3, 1:3) * traslatedCylinderInertia3 * (A_b_3(1:3, 1:3))' * j_angular_link3;

B = simplify(element1 + element2 + element3); % symbolic inertia matrix
kinetic_energy = 0.5 * ([d_theta1; d_theta2; d_theta3])' * B * [d_theta1; d_theta2; d_theta3]; % symbolic kinetic energy

disp("Kinetic energy: ");
disp(double(subs(kinetic_energy, [theta1 theta2 theta3 d_theta1 d_theta2 d_theta3], [jointValues, velocityValues]))); % numeric kinetic energy

%% potential energy

gravity = [0; 0; -9.81];

pot1 = mass1 * gravity.' * com1;
pot2 = mass2 * gravity.' * com2;
pot3 = mass3 * gravity.' * com3;

potential_energy =- (pot1 + pot2 + pot3); % symbolic potential energy

disp("Potential energy: ");
disp(double(subs(potential_energy, [theta1 theta2 theta3], jointValues))); % numeric potential energy

%% lagrange equation of motion with derivatives

B_time_dependent = subs(B, [theta1 theta2 theta3], [th1(t) th2(t) th3(t)]);
B_dot = diff(B_time_dependent, t);

potential_energy = subs(potential_energy, [theta1 theta2 theta3 d_theta1 d_theta2 d_theta3], [th1(t) th2(t) th3(t) d_th1(t) d_th2(t) d_th3(t)]);

first = B_time_dependent * [dd_th1(t); dd_th2(t); dd_th3(t)];
second = B_dot * [d_th1(t); d_th2(t); d_th3(t)];

product = ([d_th1(t); d_th2(t); d_th3(t)]).' * B_time_dependent * [d_th1(t); d_th2(t); d_th3(t)];

third = 0.5 * [diff(product, th1(t)), diff(product, th2(t)), diff(product, th3(t))].';
fourth = ([diff(potential_energy, th1(t)), diff(potential_energy, th2(t)), diff(potential_energy, th3(t))]).';

tau_with_derivatives = simplify(first + second - third + fourth);
tau_with_derivatives = simplify(subs(tau_with_derivatives, [diff(th1(t), t), diff(th2(t), t), diff(th3(t), t)], [d_th1(t), d_th2(t), d_th3(t)]));
tau_with_derivatives = simplify(subs(tau_with_derivatives, [th1(t) th2(t) th3(t) d_th1(t) d_th2(t) d_th3(t) dd_th1(t) dd_th2(t) dd_th3(t)], [theta1 theta2 theta3 d_theta1 d_theta2 d_theta3 dd_theta1 dd_theta2 dd_theta3]));

disp("tau with derivatives:");
disp(double(subs(tau_with_derivatives, [theta1 theta2 theta3 d_theta1 d_theta2 d_theta3 dd_theta1 dd_theta2 dd_theta3], [jointValues velocityValues accelerationValues]))); % numeric tau with derivatives

%% Lagrange equation of motion without derivatives

C_expanded = sym(zeros(3, 3, 3));
q = {theta1 theta2 theta3};
q_dot = {d_theta1 d_theta2 d_theta3};

for i = 1:3

    for j = 1:3

        for k = 1:3
            C_expanded(i, j, k) = 0.5 * (diff(B(i, j), q(k)) + diff(B(i, k), q(j)) - diff(B(j, k), q(i))) * q_dot(k);
        end

    end

end

C = simplify(sum(C_expanded, 3));

% Gravity Vector
g1 = mass1 * gravity' * j_linear_link1(1:3, 1) + mass2 * gravity' * j_linear_link2(1:3, 1) + mass3 * gravity' * j_linear_link3(1:3, 1);
g2 = mass1 * gravity' * j_linear_link1(1:3, 2) + mass2 * gravity' * j_linear_link2(1:3, 2) + mass3 * gravity' * j_linear_link3(1:3, 2);
g3 = mass1 * gravity' * j_linear_link1(1:3, 3) + mass2 * gravity' * j_linear_link2(1:3, 3) + mass3 * gravity' * j_linear_link3(1:3, 3);
g = [g1; g2; g3];

tau_without_derivatives = B * [dd_theta1; dd_theta2; dd_theta3] + C * [d_theta1; d_theta2; d_theta3] - g;

disp("tau without derivatives:");
disp(double(subs(tau_without_derivatives, [theta1 theta2 theta3 d_theta1 d_theta2 d_theta3 dd_theta1 dd_theta2 dd_theta3], [jointValues velocityValues accelerationValues]))); % numeric tau without derivatives

%% robot structure

robotStructure.DOF = 3;
robotStructure.q = [theta1 theta2 theta3];
robotStructure.dq = [d_theta1 d_theta2 d_theta3];
robotStructure.EE = A_b_ee(1:3, 4);
robotStructure.func.EE = matlabFunction(robotStructure.EE, 'Vars', robotStructure.q);
robotStructure.rotPosEE = A_b_ee;
robotStructure.func.rotPosEE = matlabFunction(robotStructure.rotPosEE, 'Vars', robotStructure.q);
robotStructure.gravity = gravity;
robotStructure.JointValues = jointValues;
robotStructure.VelocityValues = velocityValues;
robotStructure.AccelerationValues = accelerationValues;
robotStructure.DHTable = DH_matrix;

robotStructure.Jacobian.Link1 = j_link1; % symbolic
robotStructure.Jacobian.Link2 = j_link2; % symbolic
robotStructure.Jacobian.Link3 = j_link3; % symbolic

robotStructure.Jacobian.GeometricJacobian = J_geometric; % symbolic
robotStructure.func.GeometricJacobian = matlabFunction(J_geometric, 'Vars', robotStructure.q);

robotStructure.AnalitycalJacobian = J_analytical(1:3, 1:3); % symbolic
robotStructure.func.AnalitycalJacobian = matlabFunction(J_analytical(1:3, 1:3), 'Vars', robotStructure.q);
robotStructure.func.AnalitycalJacobianComplete = matlabFunction(J_analytical, 'Vars', robotStructure.q);



robotStructure.AnalitycalJacobianTranspose = J_analytical_transpose(1:3, 1:3); % symbolic
robotStructure.func.AnalitycalJacobianTranspose = matlabFunction(J_analytical_transpose(1:3, 1:3), 'Vars', robotStructure.q);
robotStructure.AnalitycalJacobianDot = J_analitycal_dot(1:3, 1:3); % symbolic
robotStructure.func.AnalitycalJacobianDot = matlabFunction(J_analitycal_dot(1:3, 1:3), 'Vars', [robotStructure.q, robotStructure.dq]);
robotStructure.AnalitycalJacobianInverse = J_analytical_inverse(1:3, 1:3); % symbolic
robotStructure.func.AnalitycalJacobianInverse = matlabFunction(J_analytical_inverse(1:3, 1:3), 'Vars', robotStructure.q);

robotStructure.KineticEnergy = kinetic_energy; % symbolic
robotStructure.PotentialEnergy = potential_energy; % symbolic

robotStructure.GMatrix = g; % symbolic
robotStructure.func.Gmatrix = matlabFunction(g, 'Vars', robotStructure.q);

robotStructure.BMatrix = B; % symbolic
robotStructure.func.Bmatrix = matlabFunction(B, 'Vars', robotStructure.q);

robotStructure.CMatrix = C; % symbolic
robotStructure.func.Cmatrix = matlabFunction(C, 'Vars', [robotStructure.q, robotStructure.dq]);

robotStructure.Inertia.Link1 = traslatedCylinderInertia1;
robotStructure.Inertia.Link2 = traslatedPrismInertia2;
robotStructure.Inertia.Link3 = traslatedCylinderInertia3;
robotStructure.CenterOfMass.Link1 = com1; % symbolic
robotStructure.CenterOfMass.Link2 = com2; % symbolic
robotStructure.CenterOfMass.Link3 = com3; % symbolic

robotStructure.RobotUrdf = robot; % from urdf file

save('Simulink/robot.mat', 'robotStructure');

%% for simulink, test for the plane

% % plane
% a = 0; b = 1; c = 0; d = 0.4;

% pos_ee = subs(A_b_ee(1:3, 4), [theta1, theta2, theta3], jointValues);

% x0 = pos_ee(1);
% y0 = pos_ee(2);
% z0 = pos_ee(3);

% % distance
% distance = (a * x0 + b * y0 + c * z0 + d) / sqrt(a ^ 2 + b ^ 2 + c ^ 2);

% disp(['distance: ', num2str(double(distance))]);

% figure(1);

% show(robot, config);
% hold on

% [X, Z] = meshgrid(-1:0.1:1, -1:0.1:1);
% Y = -d * ones(size(X));

% hold on;
% surf(X, Y, Z);

%% active compliance control

% Transformation matrix ZYZ Euler angles
T_zyz = [0, -sin(phi), cos(phi) * sin(gamma);
         0, cos(phi), sin(phi) * sin(gamma);
         1, 0, cos(gamma)];

Ta = blkdiag(eye(3), T_zyz);

% end effector frame
Te = A_b_ee;

desired_position = [0.4 -0.54 0.15]; % across the wall. robot'll reach this position with theta1=0, theta2=0, theta3=0
Td = [eye(3), desired_position'; 0 0 0 1];

% desired position in the end effector frame
T_D_e = [Td(1:3, 1:3)' * Te(1:3, 1:3), Td(1:3, 1:3)' * (Te(1:3, 4) - Td(1:3, 4));
         0 0 0 1];

Jad = inv(Ta) * blkdiag(Td(1:3, 1:3)', Td(1:3, 1:3)') * J_analytical; %#ok<*MINV>
Jad_transpose = Jad.';
