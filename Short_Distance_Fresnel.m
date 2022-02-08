% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance

% Parameters; units mm
L = 250e-3; lambda = 490e-6; k = 2*pi/lambda; R = 1e-2;
M = 1024; % samples


% Define spatial axes
dx = L/M;
x = -L/2:dx:L/2-dx;
y = x;
[X,Y] = meshgrid(x,y);

% Define frequency axes
fMax = 1/(2*dx);
df = 1/L;
fx = -fMax:df:fMax-df;
fy=fx;
[FX,FY] = meshgrid(fx,fy);

% Define initial field
% field = circularAperture(L, R, M, 0, 0);
field = exp(-4*log(2)/R^2*(X.^2+Y.^2));

%select propagation distance
w0 = fwhm2D(field, x, y);
z = 0.5*k*(w0(1)^2);
%sprintf("Rayleigh Distance Z_r = %.3e", z)
%z = 200e-2; % propagation distance

% Define Fresnel Propagtor
H = exp(-1i*pi*lambda*z*(FX.^2 + FY.^2));

% Propagate
ft = fft2(field);
proppedFt = ft .* fftshift(H);
propped = ifft2(proppedFt);

%calculate FWHM
fwhm_source = fwhm2D(abs(field).^2, x, y);
fwhm_propped = fwhm2D(abs(propped).^2, x, y);
x_ratio = fwhm_propped(1) / fwhm_source(1);
y_ratio = fwhm_propped(2) / fwhm_source(2);
sprintf("source X FWHM: %.3f\n" + ...
    "propagated X FWHM: %.3f\n" + ...
    "X FWHM ratio: %.3f\n" + ...
    "source Y FWHM: %.3f\n" + ...
    "propagated Y FWHM: %.3f\n" + ...
    "Y FWHM ratio: %.3f", [fwhm_source(1), fwhm_propped(1), ...
    x_ratio, fwhm_source(2), fwhm_propped(2), y_ratio])

% Plot
subplot(1,3,1);
imagesc(abs(field).^2);
title(sprintf('Source (FWHM=%.3f)', fwhm_source(1)));
axis('square');
colormap('gray');

subplot(1,3,2);
imagesc(real(H).*abs(fftshift(ft)));
title('Fresnel Propagator Sampling')
axis('square');
colormap('gray');

subplot(1,3,3);
imagesc(abs(propped).^2);
title(sprintf('Propagated (FWHM=%.3f)', fwhm_propped(1)));
axis('square');
colormap('gray');

function fwhm_res = fwhm2D(plane, x, y)
    %get FWHM of a 2D array along central x and y axes
    [x_Midpoint, y_Midpoint] = size(plane);
    x_dist = plane(x_Midpoint, :);
    y_dist = plane(:, y_Midpoint);
    x_fwhm = fwhm(x_dist, x);
    y_fwhm = fwhm(y_dist, y);
    fwhm_res = [x_fwhm, y_fwhm];
end

function width = fwhm(distribution, coordinates)
    %get the FWHM of an input array
    %half-max is max+min/2
    hm = (max(distribution) + min(distribution))/2;
    %get indices of the first and last half-max point
    idx1 = find((distribution >= hm), 1, 'first')
    idx2 = find(distribution >= hm, 1, 'last')
    %convert to a length based on input cooridnates
    width = coordinates(idx2) - coordinates(idx1);
end

function a = rect(x)
    a = abs(x) <= .5;
end

function a = circularAperture(L, R, M, xC, yC)
    % L source plane length (m)
    % R beam radius (m)
    % M samples
    dx = L/M;
    x = -L/2:dx:L/2-dx;
    y = x;
    [X,Y] = meshgrid(x,y);
    a = rect(.5 * ((X-xC).^2 + (Y-yC).^2) / (R^2));
end