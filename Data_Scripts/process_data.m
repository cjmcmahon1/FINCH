%Script to turn some bench images into an Fresnel hologram. Specify the
%input folder, images, and relevant plotting flags, based on level of
%detail desired. 

addpath('../MATLAB_functions/'); %include helper functions
addpath('./Data_Functions/');

%Measurements from the camera (units mm)

base_folder = '../Images/Bench_Images/5-11-22/';
%tuned crop parameters to center the images
crop1 = [1 1080 1 1440];   %5mm crop params
% crop1 = [405 1005 550 1150];   %20mm crop params
% crop1 = [450 1050 530 1130];   %75mm crop params
%open images
im1 = open_im(strcat(base_folder, 'led-500um-0deg.png'));
im2 = open_im(strcat(base_folder, 'led-500um-60deg.png'));
im3 = open_im(strcat(base_folder, 'led-500um-120deg.png'));
% convert images into data structures
figure(1);
imagesc(crop(im1, crop1));

h1 = image_data_struct(crop(im1, crop1), 0);
h2 = image_data_struct(crop(im2, crop1), 2*pi/3);
h3 = image_data_struct(crop(im3, crop1), 4*pi/3);
%calculate relevant lengths based on the cropped image size
delta_y = crop1(2) - crop1(1) + 1;
y_midpt = round(delta_y / 2);
delta_x = crop1(4) - crop1(3) + 1;
x_midpt = round(delta_x / 2); 
PARAMS  = bench_params(delta_x, delta_y);
%make complex holograms
c_hol = hol_from_data([h1 h2 h3]);
%flags to specify what plots to make
movie_label = '5um-imagIFTonly-zoom-4-27';
flag_plot_hologram = true;
flag_gen_3dhol     = true;
flag_show_focus    = false;
flag_show_movie    = true;
flag_save_movie    = false;
flag_show_PSF      = false;
%make plots of the individual images + the complex hologram
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

if flag_gen_3dhol
    figure('Name', 'Generating Movie Scans')
    z_vals = linspace(-10, 10, 200); %choose z-range to propagate
    data_hol_3d = hologram3D(c_hol, z_vals, PARAMS);
    data_hol_3d = struct('intensity', ...
                  data_hol_3d.intensity, ...
                  'x', data_hol_3d.x, ...
                  'y', data_hol_3d.y, 'z', z_vals);
    %convert to movie frames
    data_frames = hologram3D_to_frames(data_hol_3d, 'LED Through Pinhole');
    close
end

if flag_show_focus
    focus = -39.19;
    f_sz = 100; %focus size (number of pixels)
    focused = gen_hol_im(c_hol, focus, PARAMS);
%     crop_params = [70 130 70 130];
    crop_params = [y_midpt-f_sz y_midpt+f_sz x_midpt-f_sz x_midpt+f_sz];
    cropped = crop_struct(focused, crop_params);
    figure('Name', 'Focused Plots');
    label = sprintf('Focus Plane (z=%.2f)', focus);
    plot_im(cropped, label);
    hold on;
    plot([cropped.x(1) cropped.x(2*f_sz)], ...
         [cropped.y(f_sz) cropped.y(f_sz)], 'LineWidth', 1.0, ... 
          'Color', 'red');
end

if flag_show_movie 
    figure('Name', 'Z-scan Movie')
    movie(data_frames, 5, 40);
end

if flag_save_movie
    if not(isfolder('../Video'))
        mkdir('../Video')
    end
    v_noLP = VideoWriter(sprintf('../Video/%s.avi', movie_label));
    open(v_noLP);
    writeVideo(v_noLP, data_frames);
    close(v_noLP);
end

if flag_show_PSF
    %PSF Plot
    hol3d_im = squeeze(abs(data_hol_3d.intensity(y_midpt,:,:)));
    hol3d_im = hol3d_im(:,:); %crop PSF for visibility
    hol3d_ft = FT(data_hol_3d);
    hol3d_ft_im = squeeze(abs(hol3d_ft.intensity(y_midpt,:,:)));
    hol3d_ft_im = hol3d_ft_im(:,:); %crop FT
    figure('Name', 'noLP 3D PSF')
    subplot(1, 2, 1);
    imagesc(data_hol_3d.z, data_hol_3d.x(:), ...
            hol3d_im);
    colormap('gray');
    xlabel('distance from z focus (mm)');
    ylabel('x (mm)');
    title('Full 3D PSF');
    colorbar();
    subplot(1, 2, 2);
    colormap('gray');
    imagesc(hol3d_ft.fz(:), hol3d_ft.fx(:), hol3d_ft_im);
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