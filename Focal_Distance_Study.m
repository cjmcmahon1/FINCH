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
PARAMS.NA = 1e-1;        %numerical aperture

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
noise = 1.0;
midpt = round(PARAMS.Mx / 2);
z_values = linspace(0, 3, num_z_vals);
propagated = zeros(PARAMS.Mx, PARAMS.My, num_z_vals);
c_hol = zeros(PARAMS.Mx, PARAMS.My, num_z_vals);
c_hol_noisy = zeros(PARAMS.Mx, PARAMS.My, num_z_vals);

%generate an incoherent image PSF for comparison
p1 = propagate_init(0, PARAMS);
inc_psf = struct('intensity', abs(p1.field).^2, 'x', p1.x, 'y', p1.y);
figure('Name', 'Incoherent Image Plots');
subplot(1, 2, 1);
plot_im(inc_psf, 'Incoherent Imaging PSF');
inc_ft = FT(inc_psf);
ft_norm = inc_ft.intensity(midpt, midpt); %intensity normalization
inc_ft.intensity = inc_ft.intensity ./ ft_norm;
subplot(1, 2, 2);
plot(inc_ft.fx, abs(inc_ft.intensity(:,midpt)));
axis('square');
title('Incoherent Imaging MTF');
xlabel('x (mm^{-1})');

% generate signal hologram for each separation
for z_idx = 1:num_z_vals
    %for each z value, compute the initial propagated fields
    tmp_propagated = propagate_init(z_values(z_idx), PARAMS);
    propagated(:,:,z_idx) = tmp_propagated.field;
    tmp_hol = cHol(z_values(z_idx), 0, PARAMS);
    c_hol(:,:,z_idx)  = tmp_hol ./ sum(abs(tmp_hol), 'all'); %normalize
    tmp_hol_noisy = cHol(z_values(z_idx), noise, PARAMS);
    noisy_norm = sum(abs(tmp_hol_noisy), 'all');
    c_hol_noisy(:,:,z_idx)  = tmp_hol_noisy ./ noisy_norm; %normalize
end
%Here, we can note that if we put the imaging plane in the exact
%center of two point sources, the effective interference pattern is as if
%we just took the square of one of the fields. Therefore, our
%computationally efficient method for getting the interference pattern is
%to just square the single propagated field at each distance. 
int_patterns = propagated.^2;
% int_patterns = c_hol;
int_patterns = int_patterns ./ sum(abs(int_patterns), [1 2]); %normalize
ft_interference = fftshift(fft(int_patterns, [], 1), 1);
figure('Name', 'MTF: Field Squared');
legend_list = [];
for z_idx =[1 3 5]
    plt_vals = abs(ft_interference(:,midpt, z_idx));
    plot(fx, plt_vals./max(plt_vals)); %normalize
    %plot(fx, plt_vals);
    legend_list = [legend_list string(sprintf('z=%.2d', z_values(z_idx)))];
    hold on;
end
legend(legend_list);
title('Normalized MTF: Field^2');
hold off;
%Now, compute the efficiency-normalized MTF, this is just normalized by the
%zero spatial frequency component. This is a common tecnhique, but not very
%informative in analyzing FINCH, since it is a computational imaging
%technique
eff_norm_MTF = ft_interference ./ ft_interference(midpt, midpt, :);
figure('Name', 'MTF: Efficiency Normalized');
legend_list = [];
for z_idx =[1 3 5]
    plt_vals = abs(eff_norm_MTF(:,midpt, z_idx));
    plot(fx, plt_vals);
    legend_list = [legend_list string(sprintf('z=%.2d', z_values(z_idx)))];
    hold on;
end
legend(legend_list);
title('Efficiency Normalized MTF: Field^2');
hold off;
%now look at the noise. We are motivated by Heintzmann's Noise-Normalized
%MTF. We want to compare our MTF with its relative noise threshold, instead
%of arbitrarily normalizing to unit area.
c_hol_noise = abs(c_hol_noisy - c_hol);
c_hol_ft = fftshift(fft(c_hol, [], 1), 1);
c_hol_noise_ft = fftshift(fft(c_hol_noise, [], 1), 1);
figure('Name', 'MTF: Noise Normalized');
legend_list = [];
for z_idx =[1 3 5]
    plt_vals = abs(c_hol(:, midpt, z_idx));
    noise_magnitude = sum(abs(c_hol_noise(:, midpt, z_idx)), 'all');
    norm_plt_vals = plt_vals ./ noise_magnitude;
    window = abs((1:PARAMS.Mx) - midpt) <=200;
    norm_plt_vals = norm_plt_vals .* window;
    plot(fx, norm_plt_vals); %normalize
    legend_list = [legend_list string(sprintf('z=%.2d', z_values(z_idx)))];
    hold on;
end
legend(legend_list);
title('Noise-Normalized MTF');
hold off;

figure('Name', 'Noise');
legend_list = [];
for z_idx =[1, 3, 5]
    plt_vals = abs(c_hol_noise(:, midpt, z_idx));
    plot(fx, plt_vals); %normalize
    legend_list = [legend_list string(sprintf('z=%.2d', z_values(z_idx)))];
    hold on;
end
legend(legend_list);
title('Noise');
hold off;
% figure('Name', 'Movie');
% hol3d_ft = fftshift(fft2(propagated.^2));
% hol3d = struct('intensity', hol3d_ft, 'x', x, 'y', y, 'z', z_values);
% hol3d_frames = hologram3D_to_frames(hol3d, "Interference Pattern for " + ...
%     "Various Propagation Distances");
% movie(hol3d_frames);

%Sanity checks that our Fresnel propagator works correctly are in
%./Test_Scripts/

%Function Definitions are in ./MATLAB_FUNCTIONS/

function PSH_intensity = cHol(z, noise, PARAMS)
    %generate a complex hologram from an imaging plane a distance z from
    %two bandlimited point sources. Imaging plane is assumed to be
    %perfectly centered
    p1 = propagate_init(z, PARAMS);
    p2 = propagate_init(-z, PARAMS);
    PSH = complex_hologram(p1, p2, 3, noise, true);
    PSH_intensity = PSH.intensity;
end