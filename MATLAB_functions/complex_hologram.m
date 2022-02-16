function result = complex_hologram(plane, num_angles, bench_params)
    %{
    Based on Brooker(2021) equation 2.  Generate a complex-valued hologram
    from a series of (num_angles) real-valued holograms.
    %}
    arguments
        plane %interference plane we get from propagate()
        num_angles %>=3 subdivisions of 2*pi to apply differing phases
        bench_params
    end
    inc = 2*pi / num_angles; %we want (num_angles) evenly separated phases
    h_sum = zeros('like', plane.field);
    for i = 1:num_angles
        prev_angle = mod(inc*(i-1), 2*pi);
        next_angle = mod(inc*(i+1), 2*pi);
        shifted_h = shifted_hologram(plane, inc*i, bench_params, 250e-3);
        phase = exp(1i * prev_angle) - exp(1i * next_angle);
        h_sum = h_sum + shifted_h.intensity .* phase;
    end
    result = struct('intensity', h_sum, 'x', plane.x, 'y', plane.y);
end