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
%make complex holograms
c_hol_5um = hol_from_data([h1_5um h2_5um h3_5um]);
c_hol_20um = hol_from_data([h1_20um h2_20um h3_20um]);
%make plots of the individual images + the complex hologram
flag_plot_hologram = false;
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
    subplot(2, 2, 4);
    plot_im(c_hol_20um, "20um abs(Complex Hologram)");
end
flag_gen_3dhol = false;
if flag_gen_3dhol
    figure('Name', 'Generating Movie Scans')
    z_vals = linspace(-15, 15, 400); %focus seems to be ~-8mm
    %5um 3D hologram focus ~-8.56mm (frame 145)
    data_hol_3d_5um = hologram3D(c_hol_5um, z_vals, PARAMS_5um);
    %convert to movie frames
    data_frames_5um = hologram3D_to_frames(data_hol_3d_5um, '5um Pinhole');
    %20um 3D hologram
    data_hol_3d_20um = hologram3D(c_hol_20um, z_vals, PARAMS_20um);
    %convert to movie frames
    data_frames_20um = hologram3D_to_frames(data_hol_3d_20um, ...
                                            '20um Pinhole');
    close
end

flag_show_focus = true;
if flag_show_focus
    focus = -8.56;
    focused_5um = gen_hol_im(c_hol_5um, focus, PARAMS_5um);
    focused_20um = gen_hol_im(c_hol_20um, focus, PARAMS_20um);
    crop_params_5um = [70 130 70 130];
    cropped_5um = crop_struct(focused_5um, crop_params_5um);
    figure('Name', 'Focused Plot 5um');
    label_5um = sprintf('Focused 5um Pinhole (z=%.2f)', focus);
    plot_im(cropped_5um, label_5um);
    figure('Name', 'Focused Plot 20um');
    label_20um = sprintf('Focused 20um Pinhole (z=%.2f)', focus);
    plot_im(focused_20um, label_20um);
end

flag_show_movie = false;
if flag_show_movie 
    figure('Name', '5um Pinhole Z-scan Movie')
    movie(data_frames_5um, 5, 40);
    figure('Name', '20um Pinhole Z-scan Movie')
    movie(data_frames_20um, 5, 40);
end
flag_save_movie = false;
if flag_save_movie
    if not(isfolder('../Video'))
        mkdir('../Video')
    end
    v_5um = VideoWriter('../Video/5um_pinhole.avi');
    open(v_5um);
    writeVideo(v_5um, data_frames_5um);
    close(v_5um);
    v_20um = VideoWriter('../Video/20um_pinhole.avi');
    open(v_20um);
    writeVideo(v_20um, data_frames_20um);
    close(v_20um);
end
flag_show_PSF = false;
if flag_show_PSF
    %5um PSF Plot
    midpt_5um = round(delta_y_5um / 2);
    hol3d_im_5um = squeeze(abs(data_hol_3d_5um.intensity(:,midpt_5um,:)));
    hol3d_im_5um = hol3d_im_5um(60:140,:); %crop PSF for visibility
    hol3d_ft_5um = FT(data_hol_3d_5um);
    hol3d_ft_im_5um = squeeze(abs(hol3d_ft_5um.intensity(:,midpt_5um,:)));
    hol3d_ft_im_5um = hol3d_ft_im_5um(:,250:end); %crop FT
    figure('Name', '5um 3D PSF')
    subplot(1, 2, 1);
    imagesc(data_hol_3d_5um.z, data_hol_3d_5um.x(60:140), ...
            hol3d_im_5um);
    colormap('gray');
    xlabel('distance from z focus (mm)');
    ylabel('x (mm)');
    title('Full 3D PSF');
    colorbar();
    subplot(1, 2, 2);
    colormap('gray');
    imagesc(hol3d_ft_5um.fz(250:end), hol3d_ft_5um.fx(60:140), ...
            hol3d_ft_im_5um);
    xlabel('f_z (mm^{-1})');
    ylabel('f_x (mm^{-1})');
    title ('FT of Full 3D PSF');
    colorbar();
    %20um PSF Plot
    midpt_20um = round(delta_y_20um / 2);
    hol3d_im_20um = squeeze(abs(data_hol_3d_20um.intensity(:,midpt_20um,:)));
    hol3d_ft_20um = FT(data_hol_3d_20um);
    hol3d_ft_im_20um = squeeze(abs(hol3d_ft_20um.intensity(:,midpt_20um,:)));
    hol3d_ft_im_20um = hol3d_ft_im_20um(60:140,250:end); %crop FT
    figure('Name', '20um 3D PSF')
    subplot(1, 2, 1);
    imagesc(data_hol_3d_20um.z, data_hol_3d_20um.x, hol3d_im_20um);
    colormap('gray');
    xlabel('distance from z focus (mm)');
    ylabel('x (mm)');
    title('Full 3D PSF');
    colorbar();
    subplot(1, 2, 2);
    colormap('gray');
    imagesc(hol3d_ft_20um.fz(250:end), hol3d_ft_20um.fx(60:140), ...
            hol3d_ft_im_20um);
    xlabel('f_z (mm^{-1})');
    ylabel('f_x (mm^{-1})');
    title ('FT of Full 3D PSF');
    colorbar();
end

function cropped = crop(im, crop_array)
    %return a cropped image based on input array of 4 indices
    cropped = im(crop_array(1):crop_array(2), ...
                 crop_array(3):crop_array(4));
end

function cropped_struct = crop_struct(im_struct, crop_array)
    if isfield(im_struct, 'x')
        x_index = 'x';
        y_index = 'y';
    elseif isfield(im_struct, 'fx')
        x_index = 'fx';
        y_index = 'fy';
    end
    crop_int = im_struct.intensity(crop_array(1):crop_array(2), ...
                                   crop_array(3):crop_array(4));
    crop_x = im_struct.(x_index);
    crop_y = im_struct.(y_index);
    crop_x = crop_x(crop_array(1):crop_array(2));
    crop_y = crop_y(crop_array(3):crop_array(4));
    cropped_struct = struct('intensity', crop_int, ...
                            x_index, crop_x, y_index, crop_y);
end