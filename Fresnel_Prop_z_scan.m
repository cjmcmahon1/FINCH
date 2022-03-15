% Modified More_Fresnel_Prop to include a slider for Z propagation.
% Requires UIControl Toolbox
% Bug in callback function breaks plot updating, will fix later
addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.Lx = 250e-3;      %x side length of input image
PARAMS.Ly = 250e-3;      %y side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.Mx = 1024;        %x samples
PARAMS.My = 1024;        %y samples
PARAMS.NA = 0.05;        %numerical aperture

%Generate fields by Fresnel propagating constant amplitude,
%circular aperture fields two different distances z1 & z2. 
%using propagate(z, parameters)
%the Brooker papers have z1~-10mm, z2~10mm
z1 = -1; %mm
z2 = 1; %mm
z_prop = -1; %mm
p1 = propagate_init(z1, PARAMS);
p2 = propagate_init(z2, PARAMS);
%create phase shifted interference patterns
shifted1 = shifted_hologram(p1, p2, 0 * pi / 3);
shifted2 = shifted_hologram(p1, p2, 2 * pi / 3);
shifted3 = shifted_hologram(p1, p2, 4 * pi / 3);
%generate the complex-valued hologram
%This will internally generate and use the phase shifted holograms above
hol = complex_hologram(p1, p2, 3);
%fresnel propagate the complex hologram backwards/forwards
hol_plane = fresnel_prop(hol.intensity, z_prop, PARAMS);
hol_struct = struct('intensity', hol_plane, 'x', hol.x, 'y', hol.y);

%scan a range of z values and look at the reconstructed hologram

z_vals = 0:0.25:1.0;
for i =1:length(z_vals)
    subplot(1, length(z_vals), i);
    sys = gen_hol_im(hol, z_vals(i), PARAMS);
    label = sprintf('z=%3d um)', z_vals(i)*1e3);
    plot_im(sys, label);
end

function hol_struct = gen_hol_im(hol, z, bench_params)
    hol_plane = fresnel_prop(hol.intensity, z, bench_params);
    hol_struct = struct('intensity', hol_plane, 'x', hol.x, 'y', hol.y);
end