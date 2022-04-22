%Script to turn some bench images into an Fresnel hologram. Images are of
%the setup with (LP) and without (noLP) the first Linear Polarized LP1

addpath('../MATLAB_functions/'); %include helper functions
addpath('./Data_Functions/');

%Measurements from the camera (units mm)

base_folder = '../Images/Bench_Images/4-21-22/';
%tuned crop parameters to center the images
crop1 = [510 710 790 990]; %20um crop param
%open images
im1 = open_im(strcat(base_folder, '5um-0deg.png'));
im2 = open_im(strcat(base_folder, '5um-60deg.png'));
im3 = open_im(strcat(base_folder, '5um-120deg.png'));
%convert images into data structures
% figure(1);
% imagesc(crop(im1, crop1));

h1 = image_data_struct(crop(im1, crop1), 0);
h2 = image_data_struct(crop(im2, crop1), 1*pi/3);
h3 = image_data_struct(crop(im3, crop1), 2*pi/3);
%calculate relevant lengths based on the cropped image size
delta_y = crop1(2) - crop1(1) + 1;
delta_x = crop1(4) - crop1(3) + 1;
PARAMS  = bench_params(delta_x, delta_y);
%make complex holograms
c_hol = hol_from_data([h1 h2 h3]);
%make plots of the individual images + the complex hologram
flag_plot_hologram = true;
if flag_plot_hologram
    %noLP hologram plot
    hfig = figure('Name', 'LED Hologram Plot');
    pos = get(hfig,'position');
    %make plot window wider
    set(hfig,'position',pos.*[0.25 0.8 2.5 1.2]); 
    subplot(2, 2, 1);
    plot_im(h1, 'Image 1: \theta = 0');
    subplot(2, 2, 2);
    plot_im(h2, 'Image 2: \theta = 2\pi/3');
    subplot(2, 2, 3);
    plot_im(h3, 'Image 3: \theta = 4\pi/3');
    subplot(2, 2, 4);
    plot_im(c_hol, "abs(Complex Hologram)");
end
flag_gen_3dhol = true;
if flag_gen_3dhol
    figure('Name', 'Generating Movie Scans')
    z_vals = linspace(-20, 0, 200); %focus seems to be ~-8mm
    %noLP 3D hologram focus ~-8.56mm (frame 145)
    data_hol_3d = hologram3D(c_hol, z_vals, PARAMS);
    %convert to movie frames
    data_frames = hologram3D_to_frames(data_hol_3d, 'LED Through Pinhole');
    close
end

flag_show_focus = false;
if flag_show_focus
    focus = -8.56;
    focused = gen_hol_im(c_hol, focus, PARAMS);
    crop_params = [70 130 70 130];
    cropped = crop_struct(focused, crop_params);
    figure('Name', 'Focused Plots');
    label = sprintf('Focused noLP Pinhole (z=%.2f)', focus);
    plot_im(cropped, label);
end

flag_show_movie = true;
if flag_show_movie 
    figure('Name', 'noLP Pinhole Z-scan Movie')
    movie(data_frames, 5, 40);
end
flag_save_movie = false;
if flag_save_movie
    if not(isfolder('../Video'))
        mkdir('../Video')
    end
    v_noLP = VideoWriter('../Video/noLP_pinhole.avi');
    open(v_noLP);
    writeVideo(v_noLP, data_frames);
    close(v_noLP);
end
flag_show_PSF = false;
if flag_show_PSF
    %noLP PSF Plot
    midpt_noLP = round(delta_y / 2);
    hol3d_im_noLP = squeeze(abs(data_hol_3d.intensity(:,midpt_noLP,:)));
    hol3d_im_noLP = hol3d_im_noLP(60:140,:); %crop PSF for visibility
    hol3d_ft_noLP = FT(data_hol_3d);
    hol3d_ft_im_noLP = squeeze(abs(hol3d_ft_noLP.intensity(:,midpt_noLP,:)));
    hol3d_ft_im_noLP = hol3d_ft_im_noLP(:,250:end); %crop FT
    figure('Name', 'noLP 3D PSF')
    subplot(1, 2, 1);
    imagesc(data_hol_3d.z, data_hol_3d.x(60:140), ...
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