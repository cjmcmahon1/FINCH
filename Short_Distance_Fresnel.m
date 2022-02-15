% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance

addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.L = 250e-3;      %side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.M = 1024;        %samples
PARAMS.NA = 0.1;        %numerical aperture

% Define spatial axes (unused)
dx = PARAMS.L/PARAMS.M;
x = -PARAMS.L/2:dx:PARAMS.L/2-dx;
y = x;
[X,Y] = meshgrid(x,y);

% Define frequency axes (unused)
fMax = 1/(2*dx);
df = 1/PARAMS.L;
fx = -fMax:df:fMax-df;
fy=fx;
[FX,FY] = meshgrid(fx,fy);

%Generate fields by Fresnel propagating constant amplitude,
%circular aperture fields two different distances z1 & z2. 
%using propagate(z, parameters)
%the Brooker papers have z1~-10mm, z2~10mm
z1 = -1; %mm
z2 = 1; %mm
p1 = propagate(z1, PARAMS);
p2 = propagate(z2, PARAMS);
%add the two fields together
interference = struct('field', p1.field + p2.field, 'x', p1.x, 'y', p1.y);
shifted1 = shifted_hologram(interference, 0 * pi / 3, PARAMS, 250e-3);
shifted2 = shifted_hologram(interference, 2 * pi / 3, PARAMS, 250e-3);
shifted3 = shifted_hologram(interference, 4 * pi / 3, PARAMS, 250e-3);
hol = complex_hologram(interference, 3, PARAMS);
% hfig = figure;
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 3 1]); %make plot window wider
% subplot(1,3,1)
% plot_im(p1, sprintf('P1 (z=%3d um)', z1*1e3))
% subplot(1,3,2)
% plot_im(p2, sprintf('P2 (z=%3d um)', z2*1e3))
% subplot(1,3,3)
% plot_im(interference, "P1 + P2")

subplot(2, 2, 1)
plot_im(shifted1, "H1")
subplot(2, 2, 2)
plot_im(shifted2, "H2")
subplot(2, 2, 3)
plot_im(shifted3, "H3")
subplot(2, 2, 4)
plot_im(hol, "Total Hologram")
%Other sanity checks that our Fresnel propagator works correctly.

%Check that gaussian beam area doubles when we propagate
%by the Rayleigh distance.
%gaussian_beam_test()

%Check that first bessel function result (Goodman 4.4.2)
%has its first zero as predicted by the bessel function on
%Wikipedia.
%bessel_function_test()

function result = complex_hologram(plane, num_angles, bench_params)
    %Based on Brooker(2021) equation 2
    arguments
        plane %interference plane we get from propagate()
        num_angles %>=3 subdivisions of 2*pi to apply differing phases
        bench_params
    end
    inc = 2*pi / num_angles;
    h_sum = zeros('like', plane.field);
    for i = 1:num_angles
        prev_angle = mod(inc*(i-1), 2*pi);
        next_angle = mod(inc*(i+1), 2*pi);
        shifted_h = shifted_hologram(plane, inc*i, bench_params, 250e-3);
        phase = exp(1i * prev_angle) - exp(1i * next_angle);
        h_sum = h_sum + shifted_h.field .* phase;
    end
    result = struct('field', h_sum, 'x', plane.x, 'y', plane.y);
end

function result = shifted_hologram(plane, theta, bench_params, rh)
    %Based on Brooker (2021) equation 2
    arguments
        plane %interference plane we get from propagate()
        theta %artificial phase shift of the interference
        bench_params
        rh = 250e-3 %maximum radius of the hologram
    end
    P = pupil_func(rh, bench_params);
    h1 = plane.field .* exp(1i * theta);
    h2 = conj(plane.field) .* exp(-1i * theta);
    field = P .* (2 + h1 + h2);
    result = struct('field', field, 'x', plane.x, 'y', plane.y);
end

function plane = pupil_func(radius, bench_params)
    arguments
        radius %mm
        bench_params
    end
    dx = bench_params.L/bench_params.M;
    x = -bench_params.L/2:dx:bench_params.L/2-dx;
    y = x;
    [X,Y] = meshgrid(x,y);
    plane = (X.^2 + Y.^2) <= radius^2;
end

%Function Definitions
function plane_struct = propagate(zf, bench_params)
    arguments
        zf %distance from focus (mm)
        bench_params
    end
    % Define frequency axes
    dx = bench_params.L/bench_params.M;
    x = -bench_params.L/2:dx:bench_params.L/2-dx;
    y = x;
    fMax = 1/(2*dx);
    df = 1/bench_params.L;
    fx = -fMax:df:fMax-df;
    fy=fx;
    [FX,FY] = meshgrid(fx,fy);
    %Assuming we have a circular aperture illuminated by a unit-amplitude
    %plane wave, the fourier transform of the field is just a circ()
    %function with radius NA/lambda (Goodman 6.2.2):
    %This should probably be normalized in some way
    cutoff_freq = bench_params.NA / bench_params.lambda;
    fq_aperture = (FY.^2 + FX.^2) < (cutoff_freq^2);
    norm = length(x); %normalize frequency by length of spatial vector
    %not sure if the above is correct
    fq_aperture = fq_aperture * norm;
    %The Fresnel propagator is:
    H = fresnel_propagator(zf, bench_params.L, ...
                           bench_params.M, bench_params.lambda);
    %To propagate, we just multiply
    proppedFt = fftshift(fq_aperture .* H);
    plane = ifftshift(ifft2(proppedFt));
    %return struct so we can plot with correct x & y axis
    plane_struct = struct('field', plane, 'x', x, 'y', y);
end