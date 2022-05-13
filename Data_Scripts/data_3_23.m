%Script to turn some bench images into an Fresnel hologram. Images are of
%the setup with (LP) and without (noLP) the first Linear Polarized LP1

addpath('../MATLAB_functions/'); %include helper functions
addpath('./Data_Functions/');

%Measurements from the camera (units mm)

base_folder = '../Images/Bench_Images/3-23-22/';
%tuned crop parameters to center the images
crop_LP = [400 600 695 895];
crop_noLP = [400 600 685 885];
%open images
im1_noLP = open_im(strcat(base_folder, 'noLP1-5um-0deg.jpg'));
im2_noLP = open_im(strcat(base_folder, 'noLP1-5um-60deg.jpg'));
im3_noLP = open_im(strcat(base_folder, 'noLP1-5um-120deg.jpg'));
im1_LP = open_im(strcat(base_folder, 'LP1-5um-0deg.jpg'));
im2_LP = open_im(strcat(base_folder, 'LP1-5um-60deg.jpg'));
im3_LP = open_im(strcat(base_folder, 'LP1-5um-120deg.jpg'));
%convert images into data structures
% figure(1);
% subplot(1, 2, 1);
% imagesc(crop(im1_noLP, crop_noLP));
% subplot(1, 2, 2);
% imagesc(crop(im1_LP, crop_LP));
h1_noLP = image_data_struct(crop(im1_noLP, crop_noLP), 0);
h2_noLP = image_data_struct(crop(im2_noLP, crop_noLP), pi/3);
h3_noLP = image_data_struct(crop(im3_noLP, crop_noLP), 2*pi/3);
h1_LP = image_data_struct(crop(im1_LP, crop_LP), 0);
h2_LP = image_data_struct(crop(im2_LP, crop_LP), pi/3);
h3_LP = image_data_struct(crop(im3_LP, crop_LP), 2*pi/3);
%calculate relevant lengths based on the cropped image size
delta_y_noLP = crop_noLP(2) - crop_noLP(1) + 1;
delta_x_noLP = crop_noLP(4) - crop_noLP(3) + 1;
PARAMS_noLP  = bench_params(delta_x_noLP, delta_y_noLP);
delta_y_LP = crop_LP(2) - crop_LP(1) + 1;
delta_x_LP = crop_LP(4) - crop_LP(3) + 1;
PARAMS_LP  = bench_params(delta_x_LP, delta_y_LP);
%make complex holograms
c_hol_noLP = hol_from_data([h1_noLP h2_noLP h3_noLP]);
c_hol_LP = hol_from_data([h1_LP h2_LP h3_LP]);
%make plots of the individual images + the complex hologram
flag_plot_hologram = true;
if flag_plot_hologram
    %noLP hologram plot
    hfig_noLP = figure('Name', 'noLP Hologram Plot');
    pos_noLP = get(hfig_noLP,'position');
    %make plot window wider
    set(hfig_noLP,'position',pos_noLP.*[0.25 0.8 2.5 1.2]); 
    subplot(2, 2, 1);
    plot_im(h1_noLP, 'noLP Image 1: \theta = 0');
    subplot(2, 2, 2);
    plot_im(h2_noLP, 'noLP Image 2: \theta = 2\pi/3');
    subplot(2, 2, 3);
    plot_im(h3_noLP, 'noLP Image 3: \theta = 4\pi/3');
    subplot(2, 2, 4);
    plot_im(c_hol_noLP, "noLP abs(Complex Hologram)");
    %LP hologram plot
    hfig_LP = figure('Name', 'LP Hologram Plot');
    pos_LP = get(hfig_LP,'position');
    %make plot window wider
    set(hfig_LP,'position',pos_LP.*[0.25 0.8 2.5 1.2]); 
    subplot(2, 2, 1);
    plot_im(h1_LP, 'LP Image 1: \theta = 0');
    subplot(2, 2, 2);
    plot_im(h2_LP, 'LP Image 2: \theta = 2\pi/3');
    subplot(2, 2, 3);
    plot_im(h3_LP, 'LP Image 3: \theta = 4\pi/3');
    subplot(2, 2, 4);
    plot_im(c_hol_LP, "LP abs(Complex Hologram)");
end
flag_gen_3dhol = false;
if flag_gen_3dhol
    figure('Name', 'Generating Movie Scans')
    z_vals = linspace(-15, 15, 400); %focus seems to be ~-8mm
    %noLP 3D hologram focus ~-8.56mm (frame 145)
    data_hol_3d_noLP = hologram3D(c_hol_noLP, z_vals, PARAMS_noLP);
    %convert to movie frames
    data_frames_noLP = hologram3D_to_frames(data_hol_3d_noLP, 'noLP Pinhole');
    %LP 3D hologram
    data_hol_3d_LP = hologram3D(c_hol_LP, z_vals, PARAMS_LP);
    %convert to movie frames
    data_frames_LP = hologram3D_to_frames(data_hol_3d_LP, ...
                                            'LP Pinhole');
    close
end

flag_show_focus = true;
if flag_show_focus
    focus = -8.56;
    focused_noLP = gen_hol_im(c_hol_noLP, focus, PARAMS_noLP);
    focused_LP = gen_hol_im(c_hol_LP, focus, PARAMS_LP);
    crop_params_noLP = [70 130 70 130];
    cropped_noLP = crop_struct(focused_noLP, crop_params_noLP);
    figure('Name', 'Focused Plots');
    subplot(1, 2, 1);
    label_noLP = sprintf('Focused noLP Pinhole (z=%.2f)', focus);
    plot_im(cropped_noLP, label_noLP);
    subplot(1, 2, 2);
    label_LP = sprintf('Focused LP Pinhole (z=%.2f)', focus);
    plot_im(focused_LP, label_LP);
end

flag_show_movie = false;
if flag_show_movie 
    figure('Name', 'noLP Pinhole Z-scan Movie')
    movie(data_frames_noLP, 5, 40);
    figure('Name', 'LP Pinhole Z-scan Movie')
    movie(data_frames_LP, 5, 40);
end
flag_save_movie = false;
if flag_save_movie
    if not(isfolder('../Video'))
        mkdir('../Video')
    end
    v_noLP = VideoWriter('../Video/noLP_pinhole.avi');
    open(v_noLP);
    writeVideo(v_noLP, data_frames_noLP);
    close(v_noLP);
    v_LP = VideoWriter('../Video/LP_pinhole.avi');
    open(v_LP);
    writeVideo(v_LP, data_frames_LP);
    close(v_LP);
end
flag_show_PSF = false;
if flag_show_PSF
    %noLP PSF Plot
    midpt_noLP = round(delta_y_noLP / 2);
    hol3d_im_noLP = squeeze(abs(data_hol_3d_noLP.intensity(:,midpt_noLP,:)));
    hol3d_im_noLP = hol3d_im_noLP(60:140,:); %crop PSF for visibility
    hol3d_ft_noLP = FT(data_hol_3d_noLP);
    hol3d_ft_im_noLP = squeeze(abs(hol3d_ft_noLP.intensity(:,midpt_noLP,:)));
    hol3d_ft_im_noLP = hol3d_ft_im_noLP(:,250:end); %crop FT
    figure('Name', 'noLP 3D PSF')
    subplot(1, 2, 1);
    imagesc(data_hol_3d_noLP.z, data_hol_3d_noLP.x(60:140), ...
            hol3d_im_noLP);
    colormap('gray');
    xlabel('distance from z focus (mm)');
    ylabel('x (mm)');
    title('Full 3D PSF');
    colorbar();
    subplot(1, 2, 2);
    colormap('gray');
    imagesc(hol3d_ft_noLP.fz(250:end), hol3d_ft_noLP.fx(60:140), ...
            hol3d_ft_im_noLP);
    xlabel('f_z (mm^{-1})');
    ylabel('f_x (mm^{-1})');
    title ('FT of Full 3D PSF');
    colorbar();
    %LP PSF Plot
    midpt_LP = round(delta_y_LP / 2);
    hol3d_im_LP = squeeze(abs(data_hol_3d_LP.intensity(:,midpt_LP,:)));
    hol3d_ft_LP = FT(data_hol_3d_LP);
    hol3d_ft_im_LP = squeeze(abs(hol3d_ft_LP.intensity(:,midpt_LP,:)));
    hol3d_ft_im_LP = hol3d_ft_im_LP(60:140,250:end); %crop FT
    figure('Name', 'LP 3D PSF')
    subplot(1, 2, 1);
    imagesc(data_hol_3d_LP.z, data_hol_3d_LP.x, hol3d_im_LP);
    colormap('gray');
    xlabel('distance from z focus (mm)');
    ylabel('x (mm)');
    title('Full 3D PSF');
    colorbar();
    subplot(1, 2, 2);
    colormap('gray');
    imagesc(hol3d_ft_LP.fz(250:end), hol3d_ft_LP.fx(60:140), ...
            hol3d_ft_im_LP);
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