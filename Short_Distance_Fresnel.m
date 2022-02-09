% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance

% Parameters; units mm
%global L lambda k R M NA;
L = 250e-3; lambda = 490e-6; k = 2*pi/lambda; R = 5e-3;
M = 1024; % samples
NA = 0.1;

% Define spatial axes
dx = L/M;
x = -L/2:dx:L/2-dx;
y = x;
[X,Y] = meshgrid(x,y);

% Define frequency axes
fMax = 1/(2*dx);
df = 1/L;
fx = -fMax:df:fMax-df;
fy=fx;
[FX,FY] = meshgrid(fx,fy);

fq_aperture = (FY.^2 + FX.^2) < (NA/lambda)^2;
ap_ft = fftshift(fft2(fq_aperture));
%radius where bessel function has first zero
r_first_0 = 1.22/2 * lambda / NA;

%select propagation distance
%z = 200e-2; % propagation distance

% Define Fresnel Propagtor
%H = exp(-1i*pi*lambda*z*(FX.^2 + FY.^2));

%gaussian_beam_test()
%plot_im(ap_ft, 'FT Aperture')
p1 = propagate(0.1, 200e-3);
p2 = propagate(0.1, 220e-3);
interference = p1 + p2;
subplot(1,3,1)
plot_im(p1, "P1")
subplot(1,3,2)
plot_im(p2, "P2")
subplot(1,3,3)
plot_im(interference, "P1 + P2")

function H = fresnel_propagator(z, L, M, lambda)
    arguments
        z
        L = 250e-3
        M = 1024
        lambda = 490e-6
    end
    % Define Fresnel Propagtor
    %z: propagataion distance
    %global L M lambda;
    dx = L/M;
    % Define frequency axes
    fMax = 1/(2*dx);
    df = 1/L;
    fx = -fMax:df:fMax-df;
    fy=fx;
    [FX,FY] = meshgrid(fx,fy);
    H = exp(-1i*pi*lambda*z*(FX.^2 + FY.^2));
end

function plot_im(image, label)
    imagesc(abs(image).^2);
    title(label);
    axis('square');
    colormap('gray');
end

function gaussian_beam_test()
    %Run to test that when we propagate a gaussian beam by the Rayleigh
    %distance Z_r, the radius of the beam (FWHM) increases by a factor
    %of sqrt(2).
    % parameters used for the test
    L = 250e-3; lambda = 490e-6; k = 2*pi/lambda; R = 5e-2;
    M = 1024; % samples
    % Define spatial axes
    dx = L/M;
    x = -L/2:dx:L/2-dx;
    y = x;
    [X,Y] = meshgrid(x,y);

    field = exp(-4*log(2)/R^2*(X.^2+Y.^2)); % Define initial field
    w0 = sqrt(R^2 / (4*log(2))); %width of beam
    z = 0.5*k*(w0^2); %propagation distance = Rayleigh distance
    fprintf("Rayleigh Distance Z_r = %.3e\n", z)
    % Get Fresnel Propagtor
    H = fresnel_propagator(z, L, M, lambda);
    % Propagate
    ft = fft2(field);
    proppedFt = ft .* fftshift(H);
    propped = ifft2(proppedFt);
    
    %calculate FWHM
    fwhm_source = fwhm2D(abs(field), x, y);
    fwhm_propped = fwhm2D(abs(propped), x, y);
    x_ratio = fwhm_propped(1) / fwhm_source(1);
    y_ratio = fwhm_propped(2) / fwhm_source(2);
    fprintf("source X FWHM:     %.3f\n" + ...
            "propagated X FWHM: %.3f\n" + ...
            "X FWHM ratio:      %.5f\n" + ...
            "source Y FWHM:     %.3f\n" + ...
            "propagated Y FWHM: %.3f\n" + ...
            "Y FWHM ratio:      %.5f\n", ...
            [fwhm_source(1), fwhm_propped(1), ...
             x_ratio, fwhm_source(2), ...
             fwhm_propped(2), y_ratio]);
    % Plot
    subplot(1,3,1);
    imagesc(abs(field).^2);
    title(sprintf('Source (FWHM=%.3f)', fwhm_source(1)));
    axis('square');
    colormap('gray');
    
    subplot(1,3,2);
    imagesc(real(H).*abs(fftshift(ft)));
    title('Fresnel Propagator Sampling');
    axis('square');
    colormap('gray');
    
    subplot(1,3,3);
    imagesc(abs(propped).^2);
    title(sprintf('Propagated (FWHM=%.3f)', fwhm_propped(1)));
    axis('square');
    colormap('gray');
end

function fwhm_res = fwhm2D(plane, x, y)
    %get FWHM of a 2D array along central x and y axes
    midpoints = (size(plane)/2);
    x_dist = plane(midpoints(1), :);
    y_dist = plane(:, midpoints(2));
    x_fwhm = fwhm(x_dist, x);
    y_fwhm = fwhm(y_dist, y);
    fwhm_res = [x_fwhm, y_fwhm];
end

function width = fwhm(distribution, coordinates)
    %get the FWHM of an input array
    %half-max is max+min/2
    hm = (max(distribution) + min(distribution))/2;
    %get indices of the first and last half-max point
    idx1 = find((distribution >= hm), 1, 'first');
    idx2 = find((distribution >= hm), 1, 'last');
    %convert to a length based on input cooridnates
    width = coordinates(idx2) - coordinates(idx1);
end

function a = rect(x)
    a = abs(x) <= .5;
end

function a = circularAperture(L, R, M, xC, yC)
    % L source plane length (m)
    % R beam radius (m)
    % M samples
    dx = L/M;
    x = -L/2:dx:L/2-dx;
    y = x;
    [X,Y] = meshgrid(x,y);
    a = rect(.5 * ((X-xC).^2 + (Y-yC).^2) / (R^2));
end

function plane = propagate(na, zf, L, M, lambda)
    arguments
        na
        zf
        L = 250e-3
        M = 1024
        lambda = 490e-6
    end
    %na: numerical aperture
    %af: distance from focus
    % Define frequency axes
    %global L M lambda;
    dx = L/M;
    fMax = 1/(2*dx);
    df = 1/L;
    fx = -fMax:df:fMax-df;
    fy=fx;
    [FX,FY] = meshgrid(fx,fy);
    %Assuming we have a circular aperture illuminated by a unit-amplitude
    %plane wave, the fourier transform of the field is just: 
    fq_aperture = (FY.^2 + FX.^2) < (na/lambda)^2;
    %The Fresnel propagator is:
    H = fresnel_propagator(zf, L, M, lambda);
    %To propagate, we just multiply
    proppedFt = fftshift(fq_aperture .* H);
    plane = ifftshift(ifft2(proppedFt));
end