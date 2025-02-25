close all;

addpath('../utils/');

% Load the URDF model
robot = importrobot('RPR_zyx.urdf', 'DataFormat', 'row');

% home configuration
config = robot.homeConfiguration;
config(1) = pi/7;
config(2) = 0.2;
config(3) = pi/4;
q1 = config(1);
q2 = config(2);
q3 = config(3);

link1 = robot.getBody("Link2");
link2 = robot.getBody("Link3");
link3 = robot.getBody("Link4");

% link 1 cylinder
mass1 = link1.Mass;
lenght1 = 0.4;
radius1 = 0.02;
[cylinderInertia1, cylinderInertiaMatrix1, traslatedCylinderInertia1, cylinderInertiaToolbox1] = inertia_cylinder(radius1, 0, lenght1, mass1);
link1.Inertia = cylinderInertiaToolbox1;

%link 2 prism
mass2 = link2.Mass;
lenght2 = 0.3;
height2 = 0.03;
depth2 = 0.03;
[prismInertia, prismInertiaMatrix, traslatedPrismInertia, prismInertiaToolbox] = inertia_prism(lenght2, height2, depth2, mass2);
link2.Inertia = prismInertiaToolbox;

% link 3 cylinder
mass3 = link3.Mass;
lenght3 = 0.24;
radius3 = 0.02;
[cylinderInertia3, cylinderInertiaMatrix3, traslatedCylinderInertia3, cylinderInertiaToolbox3] = inertia_cylinder(radius3, 0, lenght3, mass3);
link3.Inertia = cylinderInertiaToolbox3;

velocityValues = [10 2 5];
accelerationValues = [2 25 5];
robot.Gravity = [0; 0; -9.81];
% robot.Gravity = [0; 0; 0];
inverseDynamics(robot, config, velocityValues, accelerationValues)
