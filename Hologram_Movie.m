% 3-16-22
% generate a FINCH hologram and scan various z-propagations. Then compile
% the images into a single animation.

addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.Lx = 250e-3;      %x side length of input image
PARAMS.Ly = 250e-3;      %y side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.Mx = 512;        %x samples
PARAMS.My = 512;        %y samples
PARAMS.NA = 0.1;        %numerical aperture

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
%generate the complex-valued hologram
hol = complex_hologram(p1, p2, 3);

% make a 3D hologram by Fresnel propagating various z distances
z_vals = linspace(-0.75, -0.25, 200);
hol3d = hologram3D(hol, z_vals, PARAMS);
frames = hologram3D_to_frames(hol3d, ...
                              'Simulated 2 Point Source Interference');
flag_show_movie = true;
flag_save_movie = true;
if flag_show_movie
    movie(frames, 5, 5);
end
if flag_save_movie
    v = VideoWriter('./Video/simulated_data.avi');
    open(v);
    writeVideo(v, frames);
    close(v);
end