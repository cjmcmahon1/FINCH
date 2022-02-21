function H = fresnel_propagator(z, Lx, Mx, Ly, My, lambda)
    arguments
        z % propagataion distance
        Lx = 250e-3
        Mx = 1024
        Ly = 250e-3
        My = 1024
        lambda = 490e-6
    end
    % Define Fresnel Propagtor
    dx = Lx/Mx;
    dy = Ly/My;
    % Define frequency axes
    fMax_x = 1/(2*dx);
    df_x = 1/Lx;
    fx = -fMax_x:df_x:fMax_x-df_x;
    fMax_y = 1/(2*dy);
    df_y = 1/Ly;
    fy = -fMax_y:df_y:fMax_y-df_y;
    [FX,FY] = meshgrid(fx,fy);
    H = exp(-1i*pi*lambda*z*(FX.^2 + FY.^2));
end