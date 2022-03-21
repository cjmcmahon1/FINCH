%Script to turn some bench images into an Fresnel hologram. Images taken
%are of a 5 and 20 micron pinhole

addpath('../MATLAB_functions/'); %include helper functions
addpath('./Data_Functions/');

%Measurements from the camera (units mm)

base_folder = '../Images/Bench_Images/3-18-22/';
%tuned crop parameters to center the images
crop_20um = [365 565 660 860];
crop_5um = [435 635 550 750];
%open images
im1_5um = open_im(strcat(base_folder, 'im1-5um-0deg.png'));
im2_5um = open_im(strcat(base_folder, 'im2-5um-60deg.png'));
im3_5um = open_im(strcat(base_folder, 'im3-5um-120deg.png'));
im1_20um = open_im(strcat(base_folder, 'im1-20um-0deg.png'));
im2_20um = open_im(strcat(base_folder, 'im2-20um-60deg.png'));
im3_20um = open_im(strcat(base_folder, 'im3-20um-120deg.png'));
%convert images into data structures
h1_5um = image_data_struct(crop(im1_5um, crop_5um), 0);
h2_5um = image_data_struct(crop(im2_5um, crop_5um), pi/3);
h3_5um = image_data_struct(crop(im3_5um, crop_5um), 2*pi/3);
h1_20um = image_data_struct(crop(im1_20um, crop_20um), 0);
h2_20um = image_data_struct(crop(im2_20um, crop_20um), pi/3);
h3_20um = image_data_struct(crop(im3_20um, crop_20um), 2*pi/3);
%calculate relevant lengths based on the cropped image size
delta_y_5um = crop_5um(2) - crop_5um(1) + 1;
delta_x_5um = crop_5um(4) - crop_5um(3) + 1;
PARAMS_5um  = bench_params(delta_x_5um, delta_y_5um);
delta_y_20um = crop_20um(2) - crop_20um(1) + 1;
delta_x_20um = crop_20um(4) - crop_20um(3) + 1;
PARAMS_20um  = bench_params(delta_x_20um, delta_y_20um);
%make plots of the individual images + the complex hologram
flag_plot_hologram = true;
if flag_plot_hologram
    %5um hologram plot
    hfig_5um = figure('Name', '5um Hologram Plot');
    pos_5um = get(hfig_5um,'position');
    %make plot window wider
    set(hfig_5um,'position',pos_5um.*[0.25 0.8 2.5 1.2]); 
    subplot(2, 2, 1);
    plot_im(h1_5um, '5um Image 1: \theta = 0');
    subplot(2, 2, 2);
    plot_im(h2_5um, '5um Image 2: \theta = 2\pi/3');
    subplot(2, 2, 3);
    plot_im(h3_5um, '5um Image 3: \theta = 4\pi/3');
    c_hol = hol_from_data([h1 h2 h3]);
    subplot(2, 2, 4);
    plot_im(c_hol_5um, "5um abs(Complex Hologram)");
    %20um hologram plot
    hfig_20um = figure('Name', '20um Hologram Plot');
    pos_20um = get(hfig_20um,'position');
    %make plot window wider
    set(hfig_20um,'position',pos_20um.*[0.25 0.8 2.5 1.2]); 
    subplot(2, 2, 1);
    plot_im(h1_20um, '20um Image 1: \theta = 0');
    subplot(2, 2, 2);
    plot_im(h2_20um, '20um Image 2: \theta = 2\pi/3');
    subplot(2, 2, 3);
    plot_im(h3_20um, '20um Image 3: \theta = 4\pi/3');
    c_hol = hol_from_data([h1 h2 h3]);
    subplot(2, 2, 4);
    plot_im(c_hol_20um, "20um abs(Complex Hologram)");
end
% z_vals = linspace(-15, 15, 400); %focus seems to be ~-8mm
% data_hol_3d = hologram3D(c_hol, z_vals, PARAMS);
% data_frames = hologram3D_to_frames(data_hol_3d);
% midpt = round(delta_y / 2);
% hol3d_xz_im = squeeze(abs(data_hol_3d.intensity(:,midpt,:)));
% hol3d_ft = FT(data_hol_3d);
% hol3d_ft_im = squeeze(abs(hol3d_ft.intensity(:,midpt,:)));
% subplot(1, 2, 1);
% imagesc(data_hol_3d.z, data_hol_3d.x, hol3d_xz_im);
% colormap('gray');
% xlabel('distance from z focus (mm)');
% ylabel('x (mm)');
% title('Full 3D PSF');
% colorbar();
% subplot(1, 2, 2);
% colormap('gray');
% imagesc(hol3d_ft.fz, hol3d_ft.fx, hol3d_ft_im);
% xlabel('f_z (mm^{-1})');
% ylabel('f_x (mm^{-1})');
% title ('FT of Full 3D PSF');
% colorbar();
movie(data_frames, 5, 40);

function cropped = crop(im, crop_array)
    %return a cropped image based on input array of 4 indices
    cropped = im(crop_array(1):crop_array(2), ...
                 crop_array(3):crop_array(4));
end
