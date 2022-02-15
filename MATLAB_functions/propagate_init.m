function plane_struct = propagate_init(zf, bench_params)
    %{
    Propagate a constant amplitude plane wave through a circular aperture
    a distance zf (mm). This is taken from Goodman 6.2.2 and is used as a
    way to generate two test images to interfere and generate our FINCH
    hologram. Note that we exploit the fact that we already know what the
    frequency distribution of this type of propagation should look like, so
    we start in the frequency space.
    %}
    arguments
        zf %distance from focus (mm)
        bench_params
    end
    % Define frequency axes
    dx = bench_params.L/bench_params.M;
    x = -bench_params.L/2:dx:bench_params.L/2-dx;
    y = x;
    fMax = 1/(2*dx);
    df = 1/bench_params.L;
    fx = -fMax:df:fMax-df;
    fy=fx;
    [FX,FY] = meshgrid(fx,fy);
    %Assuming we have a circular aperture illuminated by a unit-amplitude
    %plane wave, the fourier transform of the field is just a circ()
    %function with radius NA/lambda (Goodman 6.2.2):
    %This should probably be normalized in some way
    cutoff_freq = bench_params.NA / bench_params.lambda;
    fq_aperture = (FY.^2 + FX.^2) < (cutoff_freq^2);
    norm = length(x); %normalize frequency by length of spatial vector
    %not sure if the above is correct
    fq_aperture = fq_aperture * norm;
    %The Fresnel propagator is:
    H = fresnel_propagator(zf, bench_params.L, ...
                           bench_params.M, bench_params.lambda);
    %To propagate, we just multiply
    proppedFt = fftshift(fq_aperture .* H);
    plane = ifftshift(ifft2(proppedFt));
    %return struct so we can plot with correct x & y axis
    plane_struct = struct('field', plane, 'x', x, 'y', y);
end