% Chris McMahon: 5/11/22
% Given actual (approximate) measured parameters, we can generate the point
% spread hologram (PSH) for a given displacement, wavelength, and pixel
% spacing. In this script, we will convolve actual in-focus images with the
% PSH to simulate what the fresnel hologram should look like.

addpath('./MATLAB_functions/'); %include helper functions
addpath('./Data_Scripts/Data_Functions/');

%load image
im_base_folder = './Images/Bench_Images/Focused_Images/';
im = open_im(strcat(im_base_folder, 'led-500um-focused.png'));
crop_param = [1 1080 160 1240];
crop_im = crop(im, crop_param);
crop_im_hol = image_data_struct(crop_im, 0);
%load a spoke pattern as an ideal image
% crop_im = load('./Images/Bench_Images/Focused_Images/I1.mat').I1;
% crop_im_hol = image_data_struct(crop_im, 0);
%other important bench measurements
dz = 10.0; %mm
mag = 1.5; %magnification of 4f setup
NA = 0.1 / 3; %numerical aperture (scaled down by 10 - sampling)

%our 2 source points have their in focus object plane dz=1.71mm apart. In
%the imaging plane, this corresponds to dz*mag^2, due to the axial
%magnification.
z1 = -dz/2; %mm
z2 = dz/2; %mm

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
%normalize PSH to 1
PSH_norm = sum(abs(PSH.intensity), 'all');
PSH.intensity = PSH.intensity ./ PSH_norm;

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
flag_PSH_info = true;
if flag_PSH_info
    % generate lots of PSH plots to check that the sampling is sufficient
    hfig = figure('Name', 'Interference and Complex Hologram');
    pos = get(hfig,'position');
    set(hfig,'position',pos.*[0.25 0.25 2.5 1.9]); %make plot window wider
    subplot(3, 3, 1)
    p1_label = sprintf('P1 Intensity Plot (z1=%3d um)', z1*1e3);
    plot_im(p1, p1_label)
    subplot(3, 3, 2)
    p2_label = sprintf('P2 Intensity Plot (z2=%3d um)', z2*1e3);
    plot_im(p2, p2_label)
    subplot(3, 3, 3)
    plot_im(interference, 'P1 + P2 Intensity')
    subplot(3, 3, 4)
    plot_im(PSH, 'Re(Complex Hologram)', 'real')
    subplot(3, 3, 5)
    plot_im(PSH, 'Im(Complex Hologram)', 'imag')
    subplot(3, 3, 6)
    plot_im(PSH, 'Abs(Complex Hologram)', 'intensity')
    subplot(3, 3, 7)
    b_prop_label_re = sprintf('Re(Fresnel Propagated z=%3d um)', z1/2*1e3);
    plot_im(back_prop, b_prop_label_re, 'real')
    subplot(3, 3, 8)
    b_prop_label_im = sprintf('Im(Fresnel Propagated z=%3d um)', z1/2*1e3);
    plot_im(back_prop, b_prop_label_im, 'imag')
    subplot(3, 3, 9)
    b_prop_label = sprintf('Abs(Fresnel Propagated z=%3d um)', z1/2*1e3);
    plot_im(back_prop, b_prop_label, 'intensity')
end

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
plot_im(crop_im_hol, 'Focused Image');
subplot(1, 3, 2);
plot_im(conv, 'Focused Image * PSH');
subplot(1, 3, 3);
conv_bp_label = sprintf('Propagated z=%.2e um', z1/2*1e3);
plot_im(conv_bp, conv_bp_label, 'intensity');

%add shot noise (poissnrnd) to the image
noise = 1e-1;
% im_noisy = crop_im + poissrnd(noise, size(crop_im, 1), size(crop_im, 2));
% im_hol_noisy = image_data_struct(im_noisy, 0);
% conv_noisy = convolve(im_hol_noisy, PSH);
% hol_difference = abs(conv_noisy.intensity - conv.intensity);
% hol_diff_struct = struct('intensity', hol_difference, ...
%                          'x', conv_noisy.x, 'y', conv_noisy.y);

%create noisy PSH to model the noise we expect from each hologram
PSH_noisy = complex_hologram(p1, p2, 3, noise);
%normalize PSH to 1
PSH_noisy_norm = sum(abs(PSH_noisy.intensity), 'all');
PSH_noisy.intensity = PSH_noisy.intensity ./ PSH_norm;
%create phase shifted holograms for plotting
im1_noisy = convolve(crop_im_hol, PSH_noisy.images(1));
im1_noisy.angle = PSH_noisy.images(1).angle;
im2_noisy = convolve(crop_im_hol, PSH_noisy.images(2));
im2_noisy.angle = PSH_noisy.images(2).angle;
im3_noisy = convolve(crop_im_hol, PSH_noisy.images(3));
im3_noisy.angle = PSH_noisy.images(3).angle;
%convolve the in-focus image with the noisy PSH
%conv_noisy = convolve(crop_im_hol, PSH_noisy);
conv_noisy = hol_from_data([im1_noisy, im2_noisy, im3_noisy]);
%compare the noisy hologram with the noisless hologram
hol_difference = abs(conv_noisy.intensity - conv.intensity);
hol_diff_struct = struct('intensity', hol_difference, ...
                         'x', conv_noisy.x, 'y', conv_noisy.y);
%propagate the noisy hologram
conv_bp_noisy_intensity = fresnel_prop(conv_noisy.intensity, z1/2, PARAMS);
conv_bp_noisy = struct('intensity', conv_bp_noisy_intensity, ...
                   'x', crop_im_hol.x, 'y', crop_im_hol.y);
%compare the fresnel propagated noisy hologram with the noisless equivalent
bp_difference = abs(conv_bp_noisy_intensity - conv_bp_intensity);
bp_diff_struct = struct('intensity', bp_difference, ...
                        'x', conv_noisy.x, 'y', conv_noisy.y);

%Compare the noiseless and noisy PSH
figure('Name', 'PSH Comparison');
subplot(2, 3, 1);
plot_im(crop_im_hol, 'Focused Image');
subplot(2, 3, 2);
plot_im(PSH, 'Abs(PSH)', 'intensity');
subplot(2, 3, 3);
plot_im(PSH_noisy, 'Abs(Noisy PSH)', 'intensity');
subplot(2, 3, 4);
plot_im(im1_noisy, 'Noisy (\\theta=0))');
subplot(2, 3, 5);
plot_im(im2_noisy, 'Noisy (\\theta=2\pi/3)');
subplot(2, 3, 6);
plot_im(im3_noisy, 'Noisy (\\theta=4\pi/3)');
%plot the noisy hologram and the difference.
figure('Name', 'Noise Comparison');
subplot(2, 3, 1);
plot_im(conv, 'Noiseless Image');
subplot(2, 3, 2);
plot_im(conv_noisy, sprintf('\\lambda = %.2e', noise));
subplot(2, 3, 3);
plot_im(hol_diff_struct, 'Hologram Difference');
subplot(2, 3, 4);
conv_bp_label = sprintf('Propagated z=%.2e um', z1/2*1e3);
plot_im(conv_bp, conv_bp_label, 'intensity');
subplot(2, 3, 5);
conv_bp_noisy_label = sprintf('Noisy Propagated z=%.2e um', z1/2*1e3);
plot_im(conv_bp_noisy, conv_bp_noisy_label, 'intensity');
subplot(2, 3, 6);
plot_im(bp_diff_struct, 'Propagated Difference');

%Function Definitions are in ./MATLAB_functions/

function cropped = crop(im, crop_array)
    %return a cropped image based on input array of 4 indices
    cropped = im(crop_array(1):crop_array(2), ...
                 crop_array(3):crop_array(4));
end
