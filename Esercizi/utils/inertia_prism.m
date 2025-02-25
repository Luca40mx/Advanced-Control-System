function [prismInertia, prismInertiaMatrix, traslatedPrismInertia, prismInertiaToolbox] = inertia_prism(a, b, c, m)

    % w.r.t center of mass
    xx = (1/12) * m * (a ^ 2 + b ^ 2);
    xy = 0;
    xz = 0;
    yy = (1/12) * m * (a ^ 2 + c ^ 2);
    yz = 0;
    zz = (1/12) * m * (b ^ 2 + c ^ 2);

    prismInertia = [xx, yy, zz, xy, xz, yz];
    prismInertiaMatrix = [xx xy xz;
                          xy yy yz;
                          xz yz zz];

    traslatedPrismInertia = prismInertiaMatrix + m * ([-a / 2 0 0] * [-a / 2; 0; 0] * eye(3) - [-a / 2; 0; 0] * [-a / 2 0 0]);
    prismInertiaToolbox = [traslatedPrismInertia(1,1), traslatedPrismInertia(2,2), traslatedPrismInertia(3,3), traslatedPrismInertia(2,3), traslatedPrismInertia(1,3), traslatedPrismInertia(1,2)];

end
