% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance. We
% want to investigate the behavior of the system when the two focal lengths
% approach each other. We utilize some nice properties of the symmetry of
% our simulation to efficiently look at the hologram and its FT.

addpath('./MATLAB_functions/'); %include helper functions

% Parameters; units mm
PARAMS = struct;
PARAMS.Lx = 500e-3;      %x side length of input image
PARAMS.Ly = 500e-3;      %y side length of input image
PARAMS.lambda = 490e-6; %wavelength
PARAMS.Mx = 512;        %x samples
PARAMS.My = 512;        %y samples
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
num_z_vals = 10;
midpt = round(PARAMS.Mx / 2);
z_values = linspace(0, 4, num_z_vals);
propagated = zeros(PARAMS.Mx, PARAMS.My, num_z_vals);
c_hol = zeros(PARAMS.Mx, PARAMS.My, num_z_vals);
for z_idx = 1:num_z_vals
    %for each z value, compute the initial propagated fields
    tmp_propagated = propagate_init(z_values(z_idx), PARAMS);
    propagated(:,:,z_idx) = tmp_propagated.field;
    tmp_hol = cHol(z_values(z_idx), PARAMS);
    c_hol(:,:,z_idx)  = tmp_hol ./ sum(abs(tmp_hol), 'all'); %normalize
end
%Here, we can note that if we put the imaging plane in the exact
%center of two point sources, the effective interference pattern is as if
%we just took the square of one of the fields. Therefore, our
%computationally efficient method for getting the interference pattern is
%to just square the single propagated field at each distance. 
int_patterns = propagated.^2;
int_patterns = int_patterns ./ sum(abs(int_patterns), [1 2]); %normalize
ft_interference = fftshift(fft(int_patterns, [], 1), 1);
figure('Name', 'MTF');
legend_list = [];
for z_idx =[1 5 10]
    plt_vals = abs(ft_interference(:,midpt, z_idx));
    %plot(fx, plt_vals./max(plt_vals)); %normalize
    plot(fx, plt_vals); 
    legend_list(end+1) = string(sprintf('z=%.2d', z_values(z_idx)));
    hold on;
end
title('Normalized FT of Interference Pattern');
hold off;
figure('Name', 'Movie');
hol3d_ft = fftshift(fft2(propagated.^2));
hol3d = struct('intensity', hol3d_ft, 'x', x, 'y', y, 'z', z_values);
hol3d_frames = hologram3D_to_frames(hol3d, "Interference Pattern for " + ...
    "Various Propagation Distances");
movie(hol3d_frames);

%Sanity checks that our Fresnel propagator works correctly are in
%./Test_Scripts/

%Function Definitions are in ./MATLAB_FUNCTIONS/

function PSH_intensity = cHol(z, PARAMS)
    %generate a complex hologram from an imaging plane a distance z from
    %two bandlimited point sources. Imaging plane is assumed to be
    %perfectly centered
    p1 = propagate_init(z, PARAMS);
    p2 = propagate_init(-z, PARAMS);
    PSH = complex_hologram(p1, p2, 3, 0, true);
    PSH_intensity = PSH.intensity;
end