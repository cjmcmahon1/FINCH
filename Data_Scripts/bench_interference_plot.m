%Script to turn some bench images into an interference pattern. This is the
%test from 2-21-22 where I just image an illuminated pinhole when going
%through the triangle interferometer. The LP angle specified is the angle
%of the polarizer, which determines which arm of the interferometer the
%beam goes through (45 degrees for both arms). The camera is placed
%somewhere in between the two focal points of the lens. 

addpath('../MATLAB_functions/'); %include helper functions

%Measurements from the camera (units mm)

base_folder = './Images/Bench_Images/2-21-22/';

arm1 = open_im(strcat(base_folder, 'arm1-lp100deg.png'));
arm2 = open_im(strcat(base_folder, 'arm2-lp10deg.png'));
interference = open_im(strcat(base_folder, 'int-lp55deg.png'));
im1 = arm1(410:560, 390:540);
im2 = arm2(410:560, 390:540);
im3 = interference(410:560, 390:540);
comp = abs(im3 - im1/2 - im2/2);
crop_params = bench_params(150, 150);
int_plane_nobkg = struct('intensity', comp, ...
                   'x', crop_params.x, 'y', crop_params.y);
im1_plane = struct('intensity', im1, 'x', crop_params.x, ...
                   'y', crop_params.y);
im2_plane = struct('intensity', im2, 'x', crop_params.x, ...
                   'y', crop_params.y);
int_plane = struct('intensity', im3, ...
                   'x', crop_params.x, 'y', crop_params.y);
subplot(2, 2, 1)
plot_im(im1_plane, "Arm 1 (LP Theta = 0°)");
subplot(2, 2, 2)
plot_im(im2_plane, "Arm 2 (LP Theta = 90°)");
subplot(2, 2, 3)
plot_im(int_plane, "Interference (LP Theta = 45°)");
subplot(2, 2, 4);
plot_im(int_plane_nobkg, "Interference (Background Subtracted)");
colorbar()
% im1_int = im2double(rgb2gray(imread('./Images/Bench_Images/0deg.png')));
% im2_int = im2double(rgb2gray(imread('./Images/Bench_Images/120deg.png')));
% im3_int = im2double(rgb2gray(imread('./Images/Bench_Images/240deg.png')));
% im1 = struct('intensity', im1_int, 'x', X, 'y', Y);
% im2 = struct('intensity', im2_int, 'x', X, 'y', Y);
% im3 = struct('intensity', im3_int, 'x', X, 'y', Y);
% %test = shifted_hologram(im1, 2*pi/3, PARAMS, 5);
% bench_hol = bench_complex_hologram([im1 im2 im3], PARAMS);
% z_prop = 1e-3; %mm
% propped = fresnel_prop(bench_hol.intensity, z_prop, PARAMS);
% propped_im = struct('intensity', propped, 'x', ...
%                     bench_hol.x, 'y', bench_hol.y);
% subplot(2, 2, 1);
% plot_im(bench_hol, 'Bench Hologram (abs)');
% subplot(2, 2, 2);
% plot_im(propped_im, sprintf('Propagated z = %3d um', z_prop*1e3));
% subplot(2, 2, 3);
% plot_im(bench_hol, 'Bench Hologram (real)', 'real');
% subplot(2, 2, 4);
% plot_im(bench_hol, 'Bench Hologram (imag)', 'imag');


function mat = open_im(filename)
    im = imread(filename);
    mat = im2double(rgb2gray(im));
end

function PARAMS = bench_params(num_x_pixels, num_y_pixels)
    dx = 3.45e-3;
    dy = dx;
%     num_x_pixels = 1080;
%     num_y_pixels = 1440;
    L_x = dx * num_x_pixels;
    L_y = dy * num_y_pixels;
    X = -L_x/2:dx:L_x/2-dx;
    Y = -L_y/2:dy:L_y/2-dy;

    PARAMS = struct;
    PARAMS.Lx = L_x;      %x side length of input image
    PARAMS.Ly = L_y;      %y side length of input image
    PARAMS.lambda = 490e-6; %wavelength
    PARAMS.Mx = num_x_pixels;        %x samples
    PARAMS.My = num_y_pixels;        %y samples
    PARAMS.NA = 0.1;        %numerical aperture
    PARAMS.x = X;
    PARAMS.y = Y;
end
