% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance

addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
bench_parameters = struct();
bench_parameters.L = 250e-3;      %side length of input image
bemch_parameters.lambda = 490e-6; %wavelength
bench_parameters.M = 1024;        %samples
bench_parameters.NA = 0.1;        %numerical aperture
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

%run check that gaussian beam area doubles when we propagate
%by the Rayleigh distance.
%gaussian_beam_test()

%run check that first bessel function result (Goodman 4.4.2)
%has its first zero as predicted by the bessel function on
%Wikipedia.
%bessel_function_test()

%Generate fields by Fresnel propagating constant amplitude,
%circular aperture fields two different distances z1 & z2. 
%using propagate(NA, z)
p1 = propagate(0.1, 400e-3, L, M, lambda);
p2 = propagate(0.1, 415e-3, L, M, lambda);
interference = struct('field', p1.field + p2.field, 'x', p1.x, 'y', p1.y);
hfig = figure;
pos = get(hfig,'position');
set(hfig,'position',pos.*[.5 1 3 1]);
subplot(1,3,1)
plot_im(p1, "P1 (z=200um)")
subplot(1,3,2)
plot_im(p2, "P2 (z=220um)")
subplot(1,3,3)
plot_im(interference, "P1 + P2")

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

function plane_struct = propagate(na, zf, L, M, lambda)
    arguments
        na              %numerical aperture
        zf              %distance from focus (mm)
        L = 250e-3      %length scale of input image (mm)
        M = 1024        %number of samples
        lambda = 490e-6 %wavelength (mm)
    end
    % Define frequency axes
    dx = L/M;
    x = -L/2:dx:L/2-dx;
    y = x;
    fMax = 1/(2*dx);
    df = 1/L;
    fx = -fMax:df:fMax-df;
    fy=fx;
    [FX,FY] = meshgrid(fx,fy);
    %Assuming we have a circular aperture illuminated by a unit-amplitude
    %plane wave, the fourier transform of the field is just a circ()
    %function with radius NA/lambda (Goodman 6.2.2):
    %This should probably be normalized in some way
    fq_aperture = (FY.^2 + FX.^2) < (na/lambda)^2;
    %The Fresnel propagator is:
    H = fresnel_propagator(zf, L, M, lambda);
    %To propagate, we just multiply
    proppedFt = fftshift(fq_aperture .* H);
    plane = ifftshift(ifft2(proppedFt));
    %return struct so we can plot with correct x & y axis
    plane_struct = struct('field', plane, 'x', x, 'y', y);
end