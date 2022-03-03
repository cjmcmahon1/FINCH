%Check that gaussian beam area doubles when we propagate
%by the Rayleigh distance.
function gaussian_beam_test()
    %Run to test that when we propagate a gaussian beam by the Rayleigh
    %distance Z_r, the radius of the beam (FWHM) increases by a factor
    %of sqrt(2).
    % parameters used for the test
    L = 250e-3; lambda = 490e-6; k = 2*pi/lambda; R = 5e-2;
    M = 1024; % samples
    % Define spatial axes
    dx = L/M;
    x = -L/2:dx:L/2-dx;
    y = x;
    [X,Y] = meshgrid(x,y);

    field = exp(-4*log(2)/R^2*(X.^2+Y.^2)); % Define initial field
    w0 = sqrt(R^2 / (4*log(2))); %width of beam
    z = 0.5*k*(w0^2); %propagation distance = Rayleigh distance
    fprintf("Rayleigh Distance Z_r = %.3e\n", z)
    % Get Fresnel Propagtor
    H = fresnel_propagator(z, L, M, lambda);
    % Propagate
    ft = fft2(field);
    proppedFt = ft .* fftshift(H);
    propped = ifft2(proppedFt);
    
    %calculate FWHM
    fwhm_source = fwhm2D(abs(field), x, y);
    fwhm_propped = fwhm2D(abs(propped), x, y);
    x_ratio = fwhm_propped(1) / fwhm_source(1);
    y_ratio = fwhm_propped(2) / fwhm_source(2);
    fprintf("source X FWHM:     %.3f\n" + ...
            "propagated X FWHM: %.3f\n" + ...
            "X FWHM ratio:      %.5f\n" + ...
            "source Y FWHM:     %.3f\n" + ...
            "propagated Y FWHM: %.3f\n" + ...
            "Y FWHM ratio:      %.5f\n", ...
            [fwhm_source(1), fwhm_propped(1), ...
             x_ratio, fwhm_source(2), ...
             fwhm_propped(2), y_ratio]);
    % Plot
    subplot(1,3,1);
    imagesc(abs(field).^2);
    title(sprintf('Source (FWHM=%.3f)', fwhm_source(1)));
    axis('square');
    colormap('gray');
    
    subplot(1,3,2);
    imagesc(real(H).*abs(fftshift(ft)));
    title('Fresnel Propagator Sampling');
    axis('square');
    colormap('gray');
    
    subplot(1,3,3);
    imagesc(abs(propped).^2);
    title(sprintf('Propagated (FWHM=%.3f)', fwhm_propped(1)));
    axis('square');
    colormap('gray');
end

function fwhm_res = fwhm2D(plane, x, y)
    %get FWHM of a 2D array along central x and y axes
    midpoints = (size(plane)/2);
    x_dist = plane(midpoints(1), :);
    y_dist = plane(:, midpoints(2));
    x_fwhm = fwhm(x_dist, x);
    y_fwhm = fwhm(y_dist, y);
    fwhm_res = [x_fwhm, y_fwhm];
end

function width = fwhm(distribution, coordinates)
    %get the FWHM of an input array
    %half-max is max+min/2
    hm = (max(distribution) + min(distribution))/2;
    %get indices of the first and last half-max point
    idx1 = find((distribution >= hm), 1, 'first');
    idx2 = find((distribution >= hm), 1, 'last');
    %convert to a length based on input cooridnates
    width = coordinates(idx2) - coordinates(idx1);
end