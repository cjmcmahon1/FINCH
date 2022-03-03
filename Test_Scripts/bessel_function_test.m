function bessel_function_test()
    %check that we do generate a bessel function whose first zero is the
    %same as the zero given by Wikipedia
    %compare with (Goodman 4.4.2)
    % Parameters; units mm
    L = 250e-3; lambda = 490e-6;
    M = 1024; % samples
    NA = 0.01;
    
    % Define spatial axes
    dx = L/M;
    x = -L/2:dx:L/2-dx;
    
    % Define frequency axes
    fMax = 1/(2*dx);
    df = 1/L;
    fx = -fMax:df:fMax-df;
    fy=fx;
    [FX,FY] = meshgrid(fx,fy);
    
    fq_aperture = (FY.^2 + FX.^2) < (NA/lambda)^2;
    ap_ft = fftshift(fft2(fq_aperture));
    aperture_struct = struct('field', ap_ft, 'x', x, 'y', x);
    %radius where bessel function has first zero
    r_first_0 = 1.22/2 * lambda / NA;
    fprintf('First Zero of Bessel Function = %.3e\n', r_first_0);
    subplot(1,2,1);
    plot_im(aperture_struct, 'FT of Aperture');
    subplot(1,2,2);
    plot(aperture_struct.x, abs(aperture_struct.field(512,:)));
end
