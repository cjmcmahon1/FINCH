% Chris McMahon: 5/11/22
% Given actual (approximate) measured parameters, we can generate the point
% spread hologram (PSH) for a given displacement, wavelength, and pixel
% spacing. In this script, we will convolve actual in-focus images with the
% PSH to simulate what the fresnel hologram should look like.

addpath('./MATLAB_functions/'); %include helper functions
addpath('./Data_Scripts/Data_Functions/');

%load image
im_base_folder = '../Images/Bench_Images/3-15-22/';
im = open_im(strcat(im_base_folder, '5um-im1-0deg.png'));

%other important bench measurements
dz = 1.71; %mm
mag = 5; %magnification of 4f setup
NA = 25./200; %numerical aperture

%our 2 source points have their in focus object plane dz=1.71mm apart. In
%the imaging plane, this corresponds to dz*mag^2, due to the axial
%magnification.
z1 = -dz/2 * mag^2; %mm
z2 = dz/2 * mag^2; %mm

%generate relevant bench parameters given image size, NA, etc.
PARAMS = bench_params(size(im, 1), size(im, 2), NA);

%Simulate 2 defocused points which are frequency limited by the NA of the
%system. (jinc functions in real space) using propagate(z, parameters)
%the Brooker papers have z1~-10mm, z2~10mm (double check is this true?)
p1 = propagate_init(z1, PARAMS);
p2 = propagate_init(z2, PARAMS);
%add the two fields together
interference = struct('field', p1.field + p2.field, 'x', p1.x, 'y', p1.y);
%create phase shifted holograms for plotting
shifted1 = shifted_hologram(p1, p2, 0 * pi / 3);
shifted2 = shifted_hologram(p1, p2, 2 * pi / 3);
shifted3 = shifted_hologram(p1, p2, 4 * pi / 3);
%generate the complex-valued hologram
hol = complex_hologram(p1, p2, 3);
%fresnel propagate the complex hologram backwards
back_plane = fresnel_prop(hol.intensity, z_back, PARAMS);
forward_plane = fresnel_prop(hol.intensity, z_forward, PARAMS);
back_prop = struct('intensity', back_plane, 'x', hol.x, 'y', hol.y);
forward_prop = struct('intensity', forward_plane, 'x', hol.x, 'y', hol.y);

%Plot P1, P2, interference, as well as the resulting complex hologram to
%check that everything is working.
hfig = figure('Name', 'Interference and Complex Hologram');
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

%now, we want to convolve hol with the image, to get our simulated image
%hologram.


%Sanity checks that our Fresnel propagator works correctly are in
%./Test_Scripts/

%Function Definitions are in ./MATLAB_functions/
