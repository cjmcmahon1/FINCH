% More Fresnel Propagation Tests
% Same code as Short_Distance_Fresnel but with tests in response to
% comments
%{
TODO:
-Plot FT of the reconstructed hologram, look at line spread for different
propagation distances
-Check that the -1mm plot properly has both images displayed when fresnel
propagating
-Make the intensity axis consistent for different images for easier
comparison
-add images to README for reference
%}
addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.Lx = 250e-3;      %x side length of input image
PARAMS.Ly = 250e-3;      %y side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.Mx = 2048;        %x samples
PARAMS.My = 2048;        %y samples
PARAMS.NA = 0.1;        %numerical aperture

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
%create phase shifted interference patterns
shifted1 = shifted_hologram(p1, p2, 0 * pi / 3);
shifted2 = shifted_hologram(p1, p2, 2 * pi / 3);
shifted3 = shifted_hologram(p1, p2, 4 * pi / 3);
%generate the complex-valued hologram
%This will internally generate and use the phase shifted holograms above
hol = complex_hologram(p1, p2, 3);
%fresnel propagate the complex hologram backwards/forwards
back_plane = fresnel_prop(hol.intensity, z_back, PARAMS);
forward_plane = fresnel_prop(hol.intensity, z_forward, PARAMS);
back_prop = struct('intensity', back_plane, 'x', hol.x, 'y', hol.y);
forward_prop = struct('intensity', forward_plane, 'x', hol.x, 'y', hol.y);
%The backwards propagation should look like two images, z1
%propagated 0mm, and z2 propagated -2mm
back_comp_plane1 = propagate_init(0, PARAMS);
back_comp_plane2 = propagate_init(-2, PARAMS);
back_comparison = struct('field', ...
    back_comp_plane1.field + back_comp_plane2.field, ...
    'x', back_comp_plane1.x, 'y', back_comp_plane1.y);
%Forward propagation should look like z1=0, z2=2mm
fwd_comp_plane1 = propagate_init(0, PARAMS);
fwd_comp_plane2 = propagate_init(2, PARAMS);
fwd_comparison = struct('intensity', ...
    abs(fwd_comp_plane1.field^2) + abs(fwd_comp_plane2.field^2), ...
    'x', fwd_comp_plane1.x, 'y', fwd_comp_plane1.y);

%plot 3 shifted holograms
fig1 = figure('Name', 'Phase Shifted Interference Patterns');
pos = get(fig1,'position');
set(fig1,'position',pos.*[0.25 0.25 1.9 1.0]); %make plot window wider
subplot(1, 3, 1);
plot_im(shifted1, '\theta = 0');
subplot(1, 3, 2);
plot_im(shifted2, '\theta = 2\pi/3');
subplot(1, 3, 3);
plot_im(shifted3, '\theta = 4 \pi/3');
%plot comparision of the two intensities with the reconstructed hologram at
%values z = -1mm and z = +1mm
hfig = figure('Name', 'Ground Truth Field vs Propagated Hologram');
pos = get(hfig,'position');
set(hfig,'position',pos.*[0.25 0.25 2.2 1.9]); %make plot window wider
subplot(2, 2, 1)
b_prop_label = sprintf('Abs(Fresnel Propagated z=%3d um)', z_back*1e3);
plot_im(back_prop, b_prop_label, 'intensity')
subplot(2, 2, 2)
b_comp_label = sprintf('Intensity of P1 + P2 at P1');
plot_im(back_comparison, b_comp_label, 'intensity')
subplot(2, 2, 3)
f_prop_label = sprintf('Abs(Fresnel Propagated z=%3d um)', z_forward*1e3);
plot_im(forward_prop, f_prop_label, 'intensity')
subplot(2, 2, 4)
f_comp_label = sprintf('Intensity of P1 + P2 at P2');
plot_im(fwd_comparison, f_comp_label, 'intensity')