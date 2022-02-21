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
    dx = bench_params.Lx/bench_params.Mx;
    x = -bench_params.Lx/2:dx:bench_params.Lx/2-dx;
    dy = bench_params.Ly/bench_params.My;
    y = -bench_params.Ly/2:dy:bench_params.Ly/2-dy;
    fMax_x = 1/(2*dx);
    fMax_y = 1/(2*dy);
    df_x = 1/bench_params.Lx;
    fx = -fMax_x:df_x:fMax_x-df_x;
    df_y = 1/bench_params.Ly;
    fy = -fMax_y:df_y:fMax_y-df_y;
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
    H = fresnel_propagator(zf, bench_params.Lx, bench_params.Mx, ...
                           bench_params.Ly, bench_params.My, ...
                           bench_params.lambda);
    %To propagate, we just multiply
    proppedFt = fftshift(fq_aperture .* H);
    plane = ifftshift(ifft2(proppedFt));
    %return struct so we can plot with correct x & y axis
    plane_struct = struct('field', plane, 'x', x, 'y', y);
end