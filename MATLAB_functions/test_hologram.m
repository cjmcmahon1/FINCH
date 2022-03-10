% Test github code

% Simulation parameters
bench_params = struct;
bench_params.Lx = 100e-3; bench_params.Ly = 100e-3;
bench_params.Mx = 2000; bench_params.My = 2000;
bench_params.NA = .5;
bench_params.lambda = .5e-3;
zf = 15e-3;

% Expected resolution, for testing
res = 1.22*bench_params.lambda/(2*bench_params.NA);

% Test field
test = propagate_init(zf,bench_params);
field = test.field;
I = abs(field).^2;

% Test hologram
data = shifted_hologram(test, 0, bench_params, 250e-3);
Idata = data.intensity;

% make axes for region of interest
lenROI = 500;
xCen = bench_params.Mx/2; yCen = bench_params.My/2;
xROI = xCen - lenROI/2:xCen + lenROI/2 - 1;
yROI = yCen - lenROI/2:yCen + lenROI/2 - 1;

% Plot
subplot(121)
imagesc(x(xROI),y(yROI),I(xROI,yROI));
title('Single Field');
axis('square');

subplot(122)
imagesc(x(xROI),y(yROI),Idata(xROI,yROI));
title('Interference Pattern');
axis('square');

% Check - should go to zero; form 2I(1+cos(2*phi)), where phi is phase of
% single field
disp("Minimum of hologram: " + num2str(min(Idata,[],'all')));