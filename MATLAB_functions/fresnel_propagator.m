function H = fresnel_propagator(z, L, M, lambda)
    arguments
        z % propagataion distance
        L = 250e-3
        M = 1024
        lambda = 490e-6
    end
    % Define Fresnel Propagtor
    dx = L/M;
    % Define frequency axes
    fMax = 1/(2*dx);
    df = 1/L;
    fx = -fMax:df:fMax-df;
    fy=fx;
    [FX,FY] = meshgrid(fx,fy);
    H = exp(-1i*pi*lambda*z*(FX.^2 + FY.^2));
end