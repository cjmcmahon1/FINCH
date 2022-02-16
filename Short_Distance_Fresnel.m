% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance
%{
TODO:
Normalize FT in propagate() properly
%}
addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.L = 250e-3;      %side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.M = 2048;        %samples
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
z_back = -1; %mm
z_forward = 1; %mm
p1 = propagate_init(z1, PARAMS);
p2 = propagate_init(z2, PARAMS);
%add the two fields together
interference = struct('field', p1.field + p2.field, 'x', p1.x, 'y', p1.y);
%create phase shifted holograms for plotting
shifted1 = shifted_hologram(interference, 0 * pi / 3, PARAMS, 250e-3);
shifted2 = shifted_hologram(interference, 2 * pi / 3, PARAMS, 250e-3);
shifted3 = shifted_hologram(interference, 4 * pi / 3, PARAMS, 250e-3);
%generate the complex-valued hologram
hol = complex_hologram(interference, 3, PARAMS);
%fresnel propagate the complex hologram backwards
%if this is equal to z1 or z2, then we should just see a point
back_plane = fresnel_prop(hol.intensity, z_back, PARAMS);
forward_plane = fresnel_prop(hol.intensity, z_forward, PARAMS);
back_prop = struct('intensity', back_plane, 'x', hol.x, 'y', hol.y);
forward_prop = struct('intensity', forward_plane, 'x', hol.x, 'y', hol.y);

hfig = figure;
pos = get(hfig,'position');
set(hfig,'position',pos.*[0.25 0.25 2.5 1.9]); %make plot window wider
subplot(3, 3, 1)
p1_label = sprintf("P1 Intensity Plot (z1=%3d um)", z1*1e3);
plot_im(p1, p1_label)
subplot(3, 3, 2)
p2_label = sprintf("P2 Intensity Plot (z2=%3d um)", z2*1e3);
plot_im(p2, p2_label)
subplot(3, 3, 3)
plot_im(interference, "P1 + P2 Intensity")
subplot(3, 3, 4)
plot_im(hol, "Re(Complex Hologram)", 'real')
subplot(3, 3, 5)
plot_im(hol, "Im(Complex Hologram)", 'imag')
subplot(3, 3, 6)
plot_im(hol, "Abs(Complex Hologram)", 'intensity')
subplot(3, 3, 7)
b_prop_label_re = sprintf('Re(Fresnel Propagated z=%3d um)', z_back*1e3);
plot_im(back_prop, b_prop_label_re, 'real')
subplot(3, 3, 8)
b_prop_label_im = sprintf('Im(Fresnel Propagated z=%3d um)', z_back*1e3);
plot_im(back_prop, b_prop_label_im, 'imag')
subplot(3, 3, 9)
b_prop_label = sprintf('Abs(Fresnel Propagated z=%3d um)', z_back*1e3);
plot_im(back_prop, b_prop_label, 'intensity')


%Other Plots
% hfig = figure;
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 3 1]); %make plot window wider
% subplot(1,3,1)
% plot_im(p1, sprintf('P1 (z=%3d um)', z1*1e3))
% subplot(1,3,2)
% plot_im(p2, sprintf('P2 (z=%3d um)', z2*1e3))
% subplot(1,3,3)
% plot_im(interference, "P1 + P2")
% subplot(2, 3, 5)
% plot_im(shifted2, "abs(H2) (Theta = 2*pi/3)")
% subplot(2, 3, 6)
% plot_im(shifted3, "abs(H3) (Theta = 4*pi/3)")

%Other sanity checks that our Fresnel propagator works correctly.

%Check that gaussian beam area doubles when we propagate
%by the Rayleigh distance.
%gaussian_beam_test()

%Check that first bessel function result (Goodman 4.4.2)
%has its first zero as predicted by the bessel function on
%Wikipedia.
%bessel_function_test()

%Function Definitions
function result = complex_hologram(plane, num_angles, bench_params)
    %{
    Based on Brooker(2021) equation 2.  Generate a complex-valued hologram
    from a series of (num_angles) real-valued holograms.
    %}
    arguments
        plane %interference plane we get from propagate()
        num_angles %>=3 subdivisions of 2*pi to apply differing phases
        bench_params
    end
    inc = 2*pi / num_angles; %we want (num_angles) evenly separated phases
    h_sum = zeros('like', plane.field);
    for i = 1:num_angles
        prev_angle = mod(inc*(i-1), 2*pi);
        next_angle = mod(inc*(i+1), 2*pi);
        shifted_h = shifted_hologram(plane, inc*i, bench_params, 250e-3);
        phase = exp(1i * prev_angle) - exp(1i * next_angle);
        h_sum = h_sum + shifted_h.intensity .* phase;
    end
    result = struct('intensity', h_sum, 'x', plane.x, 'y', plane.y);
end

function result = shifted_hologram(plane, theta, bench_params, rh)
    %{
    Based on Brooker (2021) equation 2. Generate a real-valued hologram
    with a given phase shift theta from an input interference image i
    ntensity. We assume the input image is of the form 
    ~exp[i/z(x.^2 + y.^2)].
    %}
    arguments
        plane %interference plane we get from propagate()
        theta %artificial phase shift of the interference
        bench_params
        rh = 250e-3 %maximum radius of the hologram
    end
    P = pupil_func(rh, bench_params);
    h1 = plane.field .* exp(1i * theta);
    h2 = conj(plane.field) .* exp(-1i * theta);
    intensity = P .* (2 + h1 + h2);
    result = struct('intensity', intensity, 'x', plane.x, 'y', plane.y);
end

function plane = pupil_func(radius, bench_params)
    %pupil function in real space
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

function propped = fresnel_prop(im, zf, bench_params)
    %{
    Propagate an image (assumed to start in real space) a distance zf.
    %}
    arguments
        im %array that we want to propagate
        zf
        bench_params
    end
    H = fresnel_propagator(zf, bench_params.L, bench_params.M, ...
                           bench_params.lambda);
    % Propagate
    ft = fft2(im);
    proppedFt = ft .* fftshift(H);
    propped = ifft2(proppedFt);
end