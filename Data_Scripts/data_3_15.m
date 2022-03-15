%Script to turn some bench images into an interference pattern. This is the
%test from 2-21-22 where I just image an illuminated pinhole when going
%through the triangle interferometer. The LP angle specified is the angle
%of the polarizer, which determines which arm of the interferometer the
%beam goes through (45 degrees for both arms). The camera is placed
%somewhere in between the two focal points of the lens. 

addpath('../MATLAB_functions/'); %include helper functions

%Measurements from the camera (units mm)

base_folder = '../Images/Bench_Images/3-15-22/';

im1 = open_im(strcat(base_folder, '5um-im1-0deg.png'));
im2 = open_im(strcat(base_folder, '5um-im2-60deg.png'));
im3 = open_im(strcat(base_folder, '5um-im3-120deg.png'));
h1 = image_data_struct(im1(460:660, 820:1020), 0);
h2 = image_data_struct(im2(460:660, 820:1020), 2*pi/3);
h3 = image_data_struct(im3(460:660, 820:1020), 4*pi/3);
PARAMS = bench_params(201, 201);
% subplot(1, 3, 1)
% plot_im(h1, 'Image 1: \theta = 0');
% subplot(1, 3, 2)
% plot_im(h2, 'Image 2: \theta = 2\pi/3');
% subplot(1, 3, 3)
% plot_im(h3, 'Image 3: \theta = 4\pi/3');
c_hol = hol_from_data([h1 h2 h3]);
% plot_im(c_hol, "Complex Hologram");
hfig = figure;
pos = get(hfig,'position');
set(hfig,'position',pos.*[0.25 0.8 2.5 1.2]); %make plot window wider
%z_vals = linspace(3.8, 4.2, 4);
z_vals = [-5 -4.05 0 4.05 5];
for i = 1:length(z_vals)
    hol = gen_hol_im(c_hol, z_vals(i), PARAMS);
    subplot(1, length(z_vals), i);
    sys = gen_hol_im(hol, z_vals(i), PARAMS);
    label = sprintf('z=%3d um)', z_vals(i)*1e3);
    plot_im(sys, label);
end

function hol_struct = gen_hol_im(hol, z, bench_params)
    hol_plane = fresnel_prop(hol.intensity, z, bench_params);
    hol_struct = struct('intensity', hol_plane, 'x', hol.x, 'y', hol.y);
end

function c_hologram = hol_from_data(image_struct_list)
    num_h = length(image_struct_list);
    int_shape = size(image_struct_list(1).intensity);
    angles = zeros(1, num_h);
    intensities = zeros(num_h, int_shape(1), int_shape(2));
    for i = 1:num_h
        angles(i) = image_struct_list(i).angle;
        intensities(i,:,:) = [image_struct_list(i).intensity];
    end
    h_sum = zeros(int_shape);
    for i = 1:num_h
        prev_angle_idx = mod(i+1, num_h)+1;
        prev_angle = angles(prev_angle_idx);
        next_angle_idx = mod(i-1, num_h)+1;
        next_angle = angles(next_angle_idx);
        phase = exp(1i * prev_angle) - exp(1i * next_angle);
        h_sum = h_sum + squeeze(intensities(i,:,:)) .* phase;
    end
    c_hologram = struct('intensity', h_sum, 'x', ...
        image_struct_list(1).x, 'y', image_struct_list(1).y);
end

function mat = open_im(filename)
    im = imread(filename);
    mat = im2double(rgb2gray(im));
end

function hol = image_data_struct(image, angle)
    im_shape = size(image);
    crop_params = bench_params(im_shape(1), im_shape(2));
    hol = struct('intensity', image, 'x', crop_params.x, ...
        'y', crop_params.y, 'angle', angle);
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
