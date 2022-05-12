% Chris McMahon: 5/11/22
% Given actual (approximate) measured parameters, we can generate the point
% spread hologram (PSH) for a given displacement, wavelength, and pixel
% spacing. In this script, we will convolve actual in-focus images with the
% PSH to simulate what the fresnel hologram should look like.

addpath('./MATLAB_functions/'); %include helper functions
addpath('./Data_Scripts/Data_Functions/');

%load image
im_base_folder = './Images/Bench_Images/Focused_Images/';
im = open_im(strcat(im_base_folder, '500um-focused.png'));
crop_param = [1 1080 160 1240];
crop_im = crop(im, crop_param);
crop_im_hol = image_data_struct(crop_im, 0);

%other important bench measurements
dz = 1.71; %mm
mag = 5; %magnification of 4f setup
NA = 25./200 * 1e-1; %numerical aperture (scaled down by 10 - sampling)

%our 2 source points have their in focus object plane dz=1.71mm apart. In
%the imaging plane, this corresponds to dz*mag^2, due to the axial
%magnification.
z1 = -dz/2 * mag^2; %mm
z2 = dz/2 * mag^2; %mm

%generate relevant bench parameters given image size, NA, etc.
PARAMS = bench_params(size(crop_im, 2), size(crop_im, 1), NA);

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
PSH = complex_hologram(p1, p2, 3);

%Fresnel propagate the complex hologram backwards. We propagate backwards
%by z1/2 and forwards by z2/2 because of the phase doubling. We expect the
%image to be at z1/2, and the mirror image would be at z2/2 if it wasn't
%removed.
back_plane = fresnel_prop(PSH.intensity, z1/2, PARAMS);
forward_plane = fresnel_prop(PSH.intensity, z2/2, PARAMS);
back_prop = struct('intensity', back_plane, 'x', PSH.x, 'y', PSH.y);
forward_prop = struct('intensity', forward_plane, 'x', PSH.x, 'y', PSH.y);

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
plot_im(PSH, "Re(Complex Hologram)", 'real')
subplot(3, 3, 5)
plot_im(PSH, "Im(Complex Hologram)", 'imag')
subplot(3, 3, 6)
plot_im(PSH, "Abs(Complex Hologram)", 'intensity')
subplot(3, 3, 7)
b_prop_label_re = sprintf('Re(Fresnel Propagated z=%3d um)', z1/2*1e3);
plot_im(back_prop, b_prop_label_re, 'real')
subplot(3, 3, 8)
b_prop_label_im = sprintf('Im(Fresnel Propagated z=%3d um)', z1/2*1e3);
plot_im(back_prop, b_prop_label_im, 'imag')
subplot(3, 3, 9)
b_prop_label = sprintf('Abs(Fresnel Propagated z=%3d um)', z1/2*1e3);
plot_im(back_prop, b_prop_label, 'intensity')

%now, we want to convolve hol with the image, to get our simulated image
%hologram.
conv = convolve(crop_im_hol, PSH);
%Fresnel propagate the resulting hologram again to see the focused image
conv_bp_intensity = fresnel_prop(conv.intensity, z1/2, PARAMS);
conv_bp = struct('intensity', conv_bp_intensity, ...
                   'x', crop_im_hol.x, 'y', crop_im_hol.y);
%plot the resulting focused image and the expected hologram
figure('Name', 'Focused Image vs Hologram');
subplot(1, 3, 1);
plot_im(crop_im_hol, 'Focused 500um Pinhole');
subplot(1, 3, 2);
plot_im(conv, 'Focused Image * PSH');
subplot(1, 3, 3);
conv_bp_label = sprintf('Abs(Fresnel Propagated z=%3d um)', z1/2*1e3);
plot_im(conv_bp, conv_bp_label, 'intensity');
%Sanity checks that our Fresnel propagator works correctly are in
%./Test_Scripts/

%Function Definitions are in ./MATLAB_functions/

function cropped = crop(im, crop_array)
    %return a cropped image based on input array of 4 indices
    cropped = im(crop_array(1):crop_array(2), ...
                 crop_array(3):crop_array(4));
end
