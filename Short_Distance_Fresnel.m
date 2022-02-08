% Built off a script from 7/15/20
% Short Distance Fresnel Prop
% Start with some amplitude distribution and propagate a short distance

% Parameters; units mm
L = 250e-3; lambda = 490e-6; k = 2*pi/lambda; R = 5e-3;
M = 1024; % samples
z = 200e-3; % propagation distance

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
    
% Define Fresnel Propagtor
H = exp(-i*pi*lambda*z*(FX.^2 + FY.^2));

% Define initial field
% field = circularAperture(L, R, M, 0, 0);
field = exp(-4*log(2)/R^2*(X.^2+Y.^2));

% Propagate
ft = fft2(field);
proppedFt = ft .* fftshift(H);
propped = ifft2(proppedFt);

% Plot
subplot(1,3,1);
imagesc(abs(field).^2);
title('Source');
axis('square');
colormap('gray');

subplot(1,3,2);
imagesc(real(H).*abs(fftshift(ft)));
title('Fresnel Propagator Sampling')
axis('square');
colormap('gray');

subplot(1,3,3);
imagesc(abs(propped).^2);
title('Propagated')
axis('square');
colormap('gray');

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