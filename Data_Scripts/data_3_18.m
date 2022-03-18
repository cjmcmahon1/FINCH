%Script to turn some bench images into an interference pattern. This is the
%test from 2-21-22 where I just image an illuminated pinhole when going
%through the triangle interferometer. The LP angle specified is the angle
%of the polarizer, which determines which arm of the interferometer the
%beam goes through (45 degrees for both arms). The camera is placed
%somewhere in between the two focal points of the lens. 

addpath('../MATLAB_functions/'); %include helper functions
addpath('./Data_Functions/');

%Measurements from the camera (units mm)

base_folder = '../Images/Bench_Images/3-18-22/';
%crop = [365 565 660 860];
crop = [435 635 550 750];
im1 = open_im(strcat(base_folder, 'im1-5um-0deg.png'));
im2 = open_im(strcat(base_folder, 'im2-5um-60deg.png'));
im3 = open_im(strcat(base_folder, 'im3-5um-120deg.png'));
h1 = image_data_struct(im1(crop(1):crop(2), crop(3):crop(4)), 0);
h2 = image_data_struct(im2(crop(1):crop(2), crop(3):crop(4)), 2*pi/3);
h3 = image_data_struct(im3(crop(1):crop(2), crop(3):crop(4)), 4*pi/3);
delta_y = crop(2) - crop(1) + 1;
delta_x = crop(4) - crop(3) + 1;
PARAMS = bench_params(delta_x, delta_y);
% subplot(1, 3, 1)
% plot_im(h1, 'Image 1: \theta = 0');
% subplot(1, 3, 2)
% plot_im(h2, 'Image 2: \theta = 2\pi/3');
% subplot(1, 3, 3)
% plot_im(h3, 'Image 3: \theta = 4\pi/3');
% c_hol = hol_from_data([h1 h2 h3]);
% plot_im(c_hol, "Complex Hologram");
% hfig = figure;
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[0.25 0.8 2.5 1.2]); %make plot window wider
z_vals = linspace(-12, -4, 300); %focus seems to be ~-8mm
data_hol_3d = hologram3D(c_hol, z_vals, PARAMS);
data_frames = hologram3D_to_frames(data_hol_3d);
movie(data_frames, 5, 40);
