% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance
%{
TODO:
Normalize FT in propagate() properly
%}
addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.Lx = 250e-3;      %x side length of input image
PARAMS.Ly = 250e-3;      %y side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.Mx = 2048;        %x samples
PARAMS.My = 2048;        %y samples
PARAMS.NA = 0.1;        %numerical aperture

% Define spatial axes (unused)
dx = PARAMS.Lx/PARAMS.Mx;
x = -PARAMS.Lx/2:dx:PARAMS.Lx/2-dx;
dy = PARAMS.Ly/PARAMS.My;
y = -PARAMS.Ly/2:dy:PARAMS.Ly/2-dy;
[X,Y] = meshgrid(x,y);

% Define frequency axes (unused)
fMax = 1/(2*dx);
df = 1/PARAMS.Lx;
fx = -fMax:df:fMax-df;
fy=fx;
[FX,FY] = meshgrid(fx,fy);

%Generate fields by Fresnel propagating constant amplitude,
%circular aperture fields two different distances z1 & z2. 
%using propagate(z, parameters)
%the Brooker papers have z1~-10mm, z2~10mm
z1 = -1; %mm
z2 = 1; %mm
z_back = -1; %mm
z_forward = 1; %mm
p1 = propagate_init(z1, PARAMS);
p2 = propagate_init(z2, PARAMS);
%add the two fields together
interference = struct('field', p1.field + p2.field, 'x', p1.x, 'y', p1.y);
%create phase shifted holograms for plotting
shifted1 = shifted_hologram(p1, p2, 0 * pi / 3);
shifted2 = shifted_hologram(p1, p2, 2 * pi / 3);
shifted3 = shifted_hologram(p1, p2, 4 * pi / 3);
%generate the complex-valued hologram
hol = complex_hologram(p1, p2, 3);
%fresnel propagate the complex hologram backwards
%if this is equal to z1 or z2, then we should just see a point
back_plane = fresnel_prop(hol.intensity, z_back, PARAMS);
forward_plane = fresnel_prop(hol.intensity, z_forward, PARAMS);
back_prop = struct('intensity', back_plane, 'x', hol.x, 'y', hol.y);
forward_prop = struct('intensity', forward_plane, 'x', hol.x, 'y', hol.y);

% hfig = figure;
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[0.25 0.25 2.5 1.9]); %make plot window wider
% subplot(3, 3, 1)
% p1_label = sprintf("P1 Intensity Plot (z1=%3d um)", z1*1e3);
% plot_im(p1, p1_label)
% subplot(3, 3, 2)
% p2_label = sprintf("P2 Intensity Plot (z2=%3d um)", z2*1e3);
% plot_im(p2, p2_label)
% subplot(3, 3, 3)
% plot_im(interference, "P1 + P2 Intensity")
% subplot(3, 3, 4)
% plot_im(hol, "Re(Complex Hologram)", 'real')
% subplot(3, 3, 5)
% plot_im(hol, "Im(Complex Hologram)", 'imag')
% subplot(3, 3, 6)
% plot_im(hol, "Abs(Complex Hologram)", 'intensity')
% subplot(3, 3, 7)
% b_prop_label_re = sprintf('Re(Fresnel Propagated z=%3d um)', z_back*1e3);
% plot_im(back_prop, b_prop_label_re, 'real')
% subplot(3, 3, 8)
% b_prop_label_im = sprintf('Im(Fresnel Propagated z=%3d um)', z_back*1e3);
% plot_im(back_prop, b_prop_label_im, 'imag')
% subplot(3, 3, 9)
% b_prop_label = sprintf('Abs(Fresnel Propagated z=%3d um)', z_back*1e3);
% plot_im(back_prop, b_prop_label, 'intensity')

% hol3d_ft_xz = ft(ft(hol3d.intensity, 1), 3);
% hol3d_ft_struct = struct('intensity', hol3d_ft_xz(:,1024,:), ...
%     'x', hol3d.x, 'y', hol3d.y);
% plot_im(hol3d_ft_struct, '3D PSF');

% plot_im(FT(back_prop), "FT of back\_prop")

% make a 3D hologram by Fresnel propagating various z distances
z_vals = linspace(-1, 0, 30);
hol3d = hologram3D(hol, z_vals, PARAMS);
frames = hologram3D_to_frames(hol3d);
movie(frames, 5, 5);

function F = hologram3D_to_frames(hologram_struct)
    z_vals = hologram_struct.z;
    hol_max = max(abs(hologram_struct.intensity(:)));
    F(length(z_vals)) = struct('cdata',[],'colormap',[]);
    for i = 1:length(z_vals)
        X = abs(squeeze(hologram_struct.intensity(:,:,i)));
%         tmp_label = sprintf('z = %3d \mum', z_vals(i)*1e3);
        imagesc(hologram_struct.x, hologram_struct.y, X./hol_max);
        colormap('gray');
        caxis([0 1]);
        axis('square');
        xlabel("x (mm)");
        ylabel("y (mm)");
        colorbar();
        drawnow
        F(i) = getframe;
        %F(i) = im2frame(x_approx, newmap);
    end

end

function hol3d_struct = hologram3D(hologram, z_vals, bench_params)
    %generate a rank-3 complex hologram by fresnel propagating along
    %z_vals (assumed to be evenly spaced)
    num_z = length(z_vals);
    hol_shape = size(hologram.intensity);
    hol3d = zeros(hol_shape(1), hol_shape(2),num_z);
    for i = 1:num_z
        hol3d(:,:,i) = fresnel_prop(hologram.intensity, ...
                                    z_vals(i), bench_params);
    end
    hol3d_struct = struct('intensity', hol3d, 'x', hologram.x, 'y', ...
                  hologram.y, 'z', z_vals);
end

function plane_struct = FT(image_struct)
    if isfield(image_struct, 'intensity')
        field_type = 'intensity';
    elseif isfield(image_struct, 'field')
        field_type = 'field';
    else
        fprintf("Struct did not have an 'intensity' or 'field' field");
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
       plane_struct.z = fz;
    end
end

%Other Plots
% hfig = figure;
% pos = get(hfig,'position');
% set(hfig,'position',pos.*[.5 1 3 1]); %make plot window wider
% subplot(1,3,1)
% plot_im(p1, sprintf('P1 (z=%3d um)', z1*1e3))
% subplot(1,3,2)
% plot_im(p2, sprintf('P2 (z=%3d um)', z2*1e3))
% subplot(1,3,3)
% plot_im(interference, "P1 + P2")
% subplot(2, 3, 5)
% plot_im(shifted2, "abs(H2) (Theta = 2*pi/3)")
% subplot(2, 3, 6)
% plot_im(shifted3, "abs(H3) (Theta = 4*pi/3)")

%Sanity checks that our Fresnel propagator works correctly are in
%./Test_Scripts/

%Function Definitions are in ./MATLAB_FUNCTIONS
%Not yet working.
function plane = zoom(image_struct)
    if isfield(image_struct, 'intensity')
        field_type = 'intensity';
    elseif isfield(image_struct, 'field')
        field_type = 'field';
    else
        fprintf("Struct did not have an 'intensity' or 'field' field");
    end
    midpt = idivide(size(image_struct.(field_type)), 2);
    buffer = idivide(size(image_struct.(field_type)), 50);
    scan = image_struct.(field_type)(midpt);
    threshold = max(scan) / 25;
    first = find(scan > threshold, 1, "first") - buffer;
    last = find(scan > threshold, 1, "last") + buffer;
    new_field = image_struct.(field_type)(first:last, first:last);
    new_x = image_struct.x(first:last);
    new_y = new_x;
    plane = struct(field_type, new_field, 'x', new_x, 'y', new_y);
end

function b=ft(a, dim)
    % 2D Fourier transform
    b=fftshift(fft(ifftshift(a, dim), dim), dim);
end

function b=ift(a, dim)
    % 2D inversion Fourier transform
    b=ifftshift(ifft(fftshift(a, dim), dim), dim);
end