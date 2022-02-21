%Script to turn some bench images into a hologram

addpath('./MATLAB_functions/'); %include helper functions

%Measurements from the camera (units mm)
dx = 3.45e-3;
dy = dx;
num_x_pixels = 1080;
num_y_pixels = 1440;
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

im1_int = im2double(rgb2gray(imread('./Images/Bench_Images/0deg.png')));
im2_int = im2double(rgb2gray(imread('./Images/Bench_Images/120deg.png')));
im3_int = im2double(rgb2gray(imread('./Images/Bench_Images/240deg.png')));
im1 = struct('intensity', im1_int, 'x', X, 'y', Y);
im2 = struct('intensity', im2_int, 'x', X, 'y', Y);
im3 = struct('intensity', im3_int, 'x', X, 'y', Y);
%test = shifted_hologram(im1, 2*pi/3, PARAMS, 5);
bench_hol = bench_complex_hologram([im1 im2 im3], PARAMS);
z_prop = 1e-3; %mm
propped = fresnel_prop(bench_hol.intensity, z_prop, PARAMS);
propped_im = struct('intensity', propped, 'x', bench_hol.x, 'y', bench_hol.y);
subplot(2, 2, 1);
plot_im(bench_hol, 'Bench Hologram (abs)');
subplot(2, 2, 2);
plot_im(propped_im, sprintf('Propagated z = %3d um', z_prop*1e3));
subplot(2, 2, 3);
plot_im(bench_hol, 'Bench Hologram (real)', 'real');
subplot(2, 2, 4);
plot_im(bench_hol, 'Bench Hologram (imag)', 'imag');