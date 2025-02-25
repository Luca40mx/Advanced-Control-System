function [cylinderInertia, cylinderInertiaMatrix, traslatedCylinderInertia, cylinderInertiaToolbox] = inertia_cylinder(a, b, h, m)

    % w.r.t center of mass
    % h: length,  a: outer radius, b: inner radius
    xx = 0.5 * m * (a ^ 2 + b ^ 2);
    xy = 0;
    xz = 0;
    yy = 0.5 * m * (3 * (a ^ 2 + b ^ 2)^2 + h ^ 2);
    yz = 0;
    zz = 0.5 * m * (3 * (a ^ 2 + b ^ 2)^2 + h ^ 2);

    cylinderInertia = [xx, yy, zz, yz, xz, xy];
    cylinderInertiaMatrix = [xx xy xz;
                             xy yy yz;
                             xz yz zz];

    traslatedCylinderInertia = cylinderInertiaMatrix + m * ([-h / 2, 0, 0] * [-h / 2; 0; 0] * eye(3) - [-h / 2; 0; 0] * [-h / 2, 0, 0]);

    cylinderInertiaToolbox = [traslatedCylinderInertia(1,1), traslatedCylinderInertia(2,2), traslatedCylinderInertia(3,3), traslatedCylinderInertia(2,3), traslatedCylinderInertia(1,3), traslatedCylinderInertia(1,2)];
end
