function propped = fresnel_prop(im, zf, bench_params)
    %{
    Propagate an image (assumed to start in real space) a distance zf.
    %}
    arguments
        im %array that we want to propagate
        zf
        bench_params
    end
    H = fresnel_propagator(zf, bench_params.Lx, bench_params.Mx, ...
                            bench_params.Ly, bench_params.My, ...
                            bench_params.lambda);
    % Propagate
    ft = fft2(im);
    proppedFt = ft .* fftshift(H);
    propped = ifft2(proppedFt);
end