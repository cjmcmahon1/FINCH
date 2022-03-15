%Script to turn some bench images into an interference pattern. This is the
%test from 2-21-22 where I just image an illuminated pinhole when going
%through the triangle interferometer. The LP angle specified is the angle
%of the polarizer, which determines which arm of the interferometer the
%beam goes through (45 degrees for both arms). The camera is placed
%somewhere in between the two focal points of the lens. 

addpath('../MATLAB_functions/'); %include helper functions
addpath('./Data_Functions/');

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