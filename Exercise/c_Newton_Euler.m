clearvars; close all;
addpath("utils \");

%%%%%%%%%%% parameters %%%%%%%%%%%
lb = 0.15;
l1 = 0.4;
l2 = 0.3;
l3 = 0.24;
L = [l1 l2 l3];
gravity = [0; 0; -9.81];
N = 3; % number of links

jointValues = [pi / 2 -0.2 pi / 3];
velocityValues = [5 2 8];
accelerationValues = [4 1 2];
fv = 0;
fc = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

robot = importrobot("RPR_zyx_LucaPonti.urdf", "urdf", DataFormat = "row");

% home configuration
config = robot.homeConfiguration;

config(1) = jointValues(1);
config(2) = jointValues(2);
config(3) = jointValues(3);

%%%% inertias
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

m = [mass1; mass2; mass3];
I = {traslatedCylinderInertia1, traslatedPrismInertia2, traslatedCylinderInertia3};

%%%%%%%%%%% Direct Kinematics %%%%%%%%%%%

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

%%%%%%%%%%% formward equation %%%%%%%%%%%
z0 = [0; 0; 1];

prev_w = [0; 0; 0];
prev_wd = [0; 0; 0];
prev_pdd = [0; 0; 0] - gravity;
% prev_pdd = [0; 0; 0];

for i = 1:N
    curr_q = jointValues(i);
    curr_vel = velocityValues(i);
    curr_acc = accelerationValues(i);

    if i == 1 % revolute joint

        R = A_0_1(1:3, 1:3);
        curr_w = R' * prev_w + R' * curr_vel * z0; % line 7
        curr_wd = R' * prev_wd + R' * (curr_acc * z0 + cross(curr_vel * prev_w, z0)); % line 8

        rici = [-L(i) / 2; 0; 0]; % Vector from origin of frame i to the center of mass of link i, negative
        rij = [L(i); 0; 0]; % Vector from origin of frame i-1 (frame of link i) to origin of i, positive

        curr_pdd = R' * prev_pdd + cross(curr_wd, rij) + cross(curr_w, cross(curr_w, rij)); % line 9

    end

    if i == 2 % prismatic joint

        R = A_1_2(1:3, 1:3);
        curr_w = R' * prev_w; % line 7
        curr_wd = R' * prev_wd; % line 8

        rici = [0; L(i) / 2; 0]; % Vector from origin of frame i to the center of mass of link i, negative
        rij = [0; -L(i); 0] - [0; jointValues(2); 0]; % Vector from origin of frame i-1 (frame of link i) to origin of i, positive

        curr_pdd = R' * prev_pdd + cross(curr_wd, rij) + cross(curr_w, cross(curr_w, rij)) + R' * curr_acc * z0 + cross(2 * curr_vel * curr_w, R' * z0); % line 9

    end

    if i == 3 % revolute joint

        R = A_2_3(1:3, 1:3);
        curr_w = R' * prev_w + R' * curr_vel * z0; % line 7
        curr_wd = R' * prev_wd + R' * (curr_acc * z0 + cross(curr_vel * prev_w, z0)); % line 8

        rici = [-L(i) / 2; 0; 0]; % Vector from origin of frame i to the center of mass of link i, negative
        rij = [L(i); 0; 0]; % Vector from origin of frame i-1 (frame of link i) to origin of i, positive

        curr_pdd = R' * prev_pdd + cross(curr_wd, rij) + cross(curr_w, cross(curr_w, rij)); % line 9

    end

    curr_pcdd = curr_pdd + cross(curr_wd, rici) + cross(curr_w, cross(curr_w, rici)); % line 10

    curr_pcdd = subs(curr_pcdd, [theta1 theta2 theta3], [jointValues(1) jointValues(2) jointValues(3)]);

    % show(robot, config); hold on; grid on; axis equal;
    % quiver3(0,0,0, curr_pcdd(1), curr_pcdd(2), curr_pcdd(3), 0.1);

    prev_w = curr_w;
    prev_wd = curr_wd;
    prev_pdd = curr_pdd;

    W{i} = curr_w; %#ok<*SAGROW>
    Wd{i} = curr_wd;
    Pdd{i} = curr_pdd;
    Pcdd{i} = curr_pcdd;
    Rij{i} = rij;
    Rici{i} = rici;
end

%%%%%%%%%%% Backward equation %%%%%%%%%%%

prev_f = [0; 0; 0];
prev_u = [0; 0; 0];

for i = N:-1:1

    curr_q = jointValues(i);
    curr_vel = velocityValues(i);
    curr_acc = accelerationValues(i);
    w = W{i};
    wd = Wd{i};
    pcdd = Pcdd{i};
    rij = Rij{i};
    rici = Rici{i};

    if i == N
        % R = eye(3);
        R = A_3_ee(1:3, 1:3);
        Rt = A_2_3(1:3, 1:3);
    end

    if i == 2
        R = A_2_3(1:3, 1:3);
        Rt = A_1_2(1:3, 1:3);
    end

    if i == 1
        R = A_1_2(1:3, 1:3);
        Rt = A_0_1(1:3, 1:3);
    end

    curr_f = R * prev_f + m(i) * pcdd; % line 7
    curr_u = cross(-curr_f, rij + rici) + R * prev_u + cross(R * prev_f, rici) + I{i} * wd + cross(w, I{i} * w); % line 8

    if i == N % revolute joint
        tau{i} = (curr_u)' * (Rt' * z0) + fv * curr_vel + fc * sign(curr_vel); % line 9
    end

    if i == 2 % prismatic joint
        tau{i} = (curr_f)' * (Rt' * z0) + fv * curr_vel + fc * sign(curr_vel); % line 9
    end

    if i == 1 % revolute joint
        tau{i} = (curr_u)' * (Rt' * z0) + fv * curr_vel + fc * sign(curr_vel); % line 9
    end

    prev_f = curr_f;
    prev_u = curr_u;

end

for i = 1:N
    disp(double(subs(tau{i}, [theta1 theta2 theta3], [jointValues(1) jointValues(2) jointValues(3)])));
end
