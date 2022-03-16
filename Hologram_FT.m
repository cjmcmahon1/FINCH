% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance

addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.Lx = 250e-3;      %x side length of input image
PARAMS.Ly = 250e-3;      %y side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.Mx = 1024;        %x samples
PARAMS.My = 1024;        %y samples
PARAMS.NA = 0.1;        %numerical aperture

%Generate fields by Fresnel propagating constant amplitude,
%circular aperture fields two different distances z1 & z2. 
%using propagate(z, parameters)
%the Brooker papers have z1~-10mm, z2~10mm
z1 = -1; %mm
z2 = 1; %mm
p1 = propagate_init(z1, PARAMS);
p2 = propagate_init(z2, PARAMS);
%generate the complex-valued hologram
hol = complex_hologram(p1, p2, 3);

% make a 3D hologram by Fresnel propagating various z distances
z_vals = linspace(-0.55, -0.45, 200);
%hol3d = hologram3D(hol, z_vals, PARAMS);
hol3d_xz_im = squeeze(abs(hol3d.intensity(:,512,:)));
subplot(1, 2, 1);
imagesc(hol3d.x, hol3d.z, hol3d_xz_im);
colormap('gray');
xlabel('x (mm)');
ylabel('z (mm)');
% hol3d_ft = FT(hol3d);
% hol3d_ft_im = squeeze(abs(hol3d_ft.intensity(:,512,:)));
% imagesc(hol3d_ft.fx, hol3d_ft.fz, hol3d_ft_im);


function plane_struct = FT(image_struct)
    if isfield(image_struct, 'intensity')
        field_type = 'intensity';
    elseif isfield(image_struct, 'field')
        field_type = 'field';
    else
        fprintf("Struct did not have an 'intensity' or 'field' field");
    end
    ft = fftshift(fftn(image_struct.(field_type)));
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