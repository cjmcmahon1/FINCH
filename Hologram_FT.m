% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance
% Update 4/5: add mirror image hologram to see all components of MTF
addpath('./MATLAB_functions/'); %include helper functions
num_pixels = 512;
midpt = num_pixels / 2;
% Parameters; units mm
PARAMS = struct;
PARAMS.Lx = 500e-3;      %x side length of input image
PARAMS.Ly = 500e-3;      %y side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.Mx = num_pixels;        %x samples
PARAMS.My = num_pixels;        %y samples
PARAMS.NA = 0.1;        %numerical aperture

%Generate fields by Fresnel propagating constant amplitude,
%circular aperture fields two different distances z1 & z2. 
%using propagate(z, parameters)
%the Brooker papers have z1~-10mm, z2~10mm
separation = 2e0; %mm
z1 = -separation/2; %mm
z2 = separation/2; %mm
p1 = propagate_init(z1, PARAMS);
p2 = propagate_init(z2, PARAMS);
const_bkg = struct('intensity', abs(p1.field).^2 + abs(p2.field).^2, ...
                   'x', p1.x, 'y', p1.y);
%generate the complex-valued hologram
hol = complex_hologram(p1, p2, 3);
hol_mirror = mirror_hologram(p1, p2, 3);
% make a 3D hologram by Fresnel propagating various z distances
%z_vals = linspace(z1/2 - 0.25, z1/2+0.25, 100);
z_vals = linspace (z1/2 - 0.25, z2/2 + 0.25, 100);
num_z_vals = size(z_vals);
num_z_vals = num_z_vals(2);

%full 3D PSF calculation
flag_calc_hologram = true;
if flag_calc_hologram
    hol3d = hologram3D(hol, z_vals, PARAMS);
    hol3d_mirror = hologram3D(hol_mirror, z_vals, PARAMS);
    hol3d_bkg = hologram3D(const_bkg, z_vals, PARAMS);
    hol3d_xz_im = squeeze(abs(hol3d.intensity(:,midpt,:)));
    hol3d_mirror_xz_im = squeeze(abs(hol3d_mirror.intensity(:,midpt,:)));
    hol3d_bkg_xz_im = squeeze(abs(hol3d_bkg.intensity(:,midpt,:)));
    hol3d_ft = FT(hol3d);
    hol3d_ft_im = squeeze(imag(hol3d_ft.intensity(:,midpt,:)));
    hol3d_mirror_ft = FT(hol3d_mirror);
    hol3d_mirror_ft_im = squeeze(imag(hol3d_mirror_ft.intensity(:,midpt,:)));
    hol3d_bkg_ft = FT(hol3d_bkg);
    hol3d_bkg_ft_im = squeeze(imag(hol3d_bkg_ft.intensity(:,midpt,:)));
end

%plot the image 3D PSF and FT
subplot(3, 2, 1);
imagesc(hol3d.z, hol3d.x, hol3d_xz_im);
colormap('gray');
xlabel('distance from z focus (mm)');
ylabel('x (mm)');
title('3D PSF');
colorbar();
subplot(3, 2, 2);
colormap('gray');
imagesc(fftshift(hol3d_ft.fz), ...
        hol3d_ft.fx, fftshift(hol3d_ft_im, 2));
xlabel('f_z (mm^{-1})');
ylabel('f_x (mm^{-1})');
title ('FT of 3D PSF');
colorbar();

%plot the mirror image 3D PSF and FT
subplot(3, 2, 3);
imagesc(hol3d_mirror.z, hol3d_mirror.x, hol3d_mirror_xz_im);
colormap('gray');
xlabel('distance from z focus (mm)');
ylabel('x (mm)');
title('Mirror Image 3D PSF');
colorbar();
subplot(3, 2, 4);
colormap('gray');
imagesc(fftshift(hol3d_mirror_ft.fz), ...
        hol3d_mirror_ft.fx, fftshift(hol3d_mirror_ft_im, 2));
xlabel('f_z (mm^{-1})');
ylabel('f_x (mm^{-1})');
title ('FT of Mirror Image 3D PSF');
colorbar();

%plot the constant background 3D PSF and FT
subplot(3, 2, 5);
imagesc(hol3d_bkg.z, hol3d_bkg.x, hol3d_bkg_xz_im);
colormap('gray');
xlabel('distance from z focus (mm)');
ylabel('x (mm)');
title('Constant Background 3D PSF');
colorbar();
subplot(3, 2, 6);
colormap('gray');
imagesc(fftshift(hol3d_bkg_ft.fz), ...
        hol3d_bkg_ft.fx, fftshift(hol3d_bkg_ft_im, 2));
xlabel('f_z (mm^{-1})');
ylabel('f_x (mm^{-1})');
title ('FT of Constant Background 3D PSF');
colorbar();

% propagate in the xz plane to speed up the calculation of 3D PSFs
% z_propped = fresnel_prop_xz(hol, midpt, z_vals, PARAMS);
% quick_ft = abs(fftshift(fft2(z_propped)));
%plot the comparison
% subplot(1, 4, 1)
% imagesc(z_vals, hol.x, abs(z_propped));
% colormap('gray');
% xlabel('z (mm)');
% ylabel('x (mm)');
% title('Quick 3D PSF');
% colorbar();
% subplot(1, 4, 2);
% imagesc(hol3d_ft.fz, hol3d_ft.fx, quick_ft);
% xlabel('f_z (mm^{-1})');
% ylabel('f_x (mm^{-1})');
% title ('FT of Quick 3D PSF');
% colorbar();

function H = fresnel_propagator_xz(z_values, yslice_idx, bench_params)
    arguments
        z_values % propagataion distance
        yslice_idx % index of desired y-slice
        bench_params % setup parameters
    end
    % Define Fresnel Propagtor
    dx = bench_params.Lx/bench_params.Mx;
    dy = bench_params.Ly/bench_params.My;
    lambda = bench_params.lambda;
    % Define frequency axes
    fMax_x = 1/(2*dx);
    df_x = 1/bench_params.Lx;
    fx = -fMax_x:df_x:fMax_x-df_x;
    fMax_y = 1/(2*dy);
    df_y = 1/bench_params.Ly;
    %get single FY value based on the y slice (this saves a lot of memory)
    FY_const = -fMax_y + (yslice_idx-1) * df_y;
    [FX,Z] = meshgrid(fx,z_values);
    quad_phase = exp(-1i*pi*lambda.*((FX.^2 + FY_const^2).*Z));
    z_phase = exp((2i*pi/lambda).*Z);
    H = quad_phase .* z_phase;
    % H = exp(2i*pi.*Z./lambda) .* exp(-1i*pi*lambda.*Z.*(FX.^2));
end

function propped = fresnel_prop_xz(hol, yslice_idx, z_values, bench_params)
    % memory-efficient function to propagate y-slice of a hologram a range
    % of distances z_values. The resulting matrix is of dimension len(x) x
    % len(z)
    arguments
        hol % 2D fresnel hologram
        yslice_idx % index of desired y-slice to propagate
        z_values % range of z values to propagate over
        bench_params
    end
    %get desired slice of hologram
    hol_slice = transpose(hol.intensity(:,yslice_idx));
    %repeat input y slice to broadcast across z values
    num_z_vals = size(z_values);
    im_xz = repmat(hol_slice, [num_z_vals(2), 1]);
    %generate fresnel propagator
    H = fresnel_propagator_xz(z_values, yslice_idx, bench_params);
    % Propagate
    im_size = size(im_xz);
    %ft only in x axis (ax=2) 
    % this is because z axis is just used for broadcasting, not propagating
    ft = fftshift(fft(im_xz, im_size(2), 2), 2);
    proppedFt = ft .* H;
    %ifft and ifftshift should also be only in x axis (ax=2)
    propped = transpose(ifft(ifftshift(proppedFt, 2), im_size(2), 2));
end

function plane_struct = FT(image_struct)
    if isfield(image_struct, 'intensity')
        field_type = 'intensity';
    elseif isfield(image_struct, 'field')
        field_type = 'field';
    else
        fprintf("Struct did not have an 'intensity' or 'field' parameter");
    end
    ft = fftshift(fftn(ifftshift(image_struct.(field_type))));
    %get correct frequency axis
    % Define spatial axes
    dx = image_struct.x(2) - image_struct.x(1);
    dy = image_struct.y(2) - image_struct.y(1);
    lx = image_struct.x(end) - image_struct.x(1);
    ly = image_struct.y(end) - image_struct.y(1);
    % Define frequency axes
    fMax_x = 1/(2*dx);
    fMax_y = 1/(2*dy);
    df_x = 1/lx;
    df_y = 1/ly;
    fx = -fMax_x:df_x:fMax_x-df_x;
    fy = -fMax_y:df_y:fMax_y-df_y;
    plane_struct = struct('intensity', ft, 'fx', fx, 'fy', fy);
    if isfield(image_struct, 'z')
       %compute z frequencies if it's a 3D hologram
       dz = image_struct.z(2) - image_struct.z(1); 
       fMax_z = 1/(2*dz);
       lz = image_struct.z(end) - image_struct.z(1);
       df_z = 1/lz;
       fz = -fMax_z:df_z:fMax_z-df_z;
       plane_struct.fz = fz;
    end
end

%Sanity checks that our Fresnel propagator works correctly are in
%./Test_Scripts/

%Function Definitions are in ./MATLAB_FUNCTIONS