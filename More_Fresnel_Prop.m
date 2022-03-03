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
%add the two fields together
interference = struct('field', p1.field + p2.field, 'x', p1.x, 'y', p1.y);
%create phase shifted holograms for plotting
shifted1 = shifted_hologram(interference, 0 * pi / 3, PARAMS, 250e-3);
shifted2 = shifted_hologram(interference, 2 * pi / 3, PARAMS, 250e-3);
shifted3 = shifted_hologram(interference, 4 * pi / 3, PARAMS, 250e-3);
%generate the complex-valued hologram
hol = complex_hologram(interference, 3, PARAMS);
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
fwd_comparison = struct('field', ...
    fwd_comp_plane1.field + fwd_comp_plane2.field, ...
    'x', fwd_comp_plane1.x, 'y', fwd_comp_plane1.y);

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