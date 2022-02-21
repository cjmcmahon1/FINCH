function result = bench_complex_hologram(plane_list, bench_params)
    %{
    Based on Brooker(2021) equation 2.  Generate a complex-valued hologram
    from a series of (num_angles) real-valued holograms.
    %}
    arguments
        plane_list %interference plane we get from propagate()
        bench_params
    end
    plane1 = plane_list(1);
    h_sum = zeros('like', plane1.intensity);
    angle_list = [0 2*pi/3 4*pi/3];
    for i = 1:length(plane_list)
        plane = plane_list(i);
        prev_angle = angle_list(mod(i-1, 3)+1);
        next_angle = angle_list(mod(i+1, 3)+1);
        angle = angle_list(mod(i, 3)+1);
        shifted_h = shifted_hologram(plane, angle, bench_params, 5);
        phase = exp(1i * prev_angle) - exp(1i * next_angle);
        h_sum = h_sum + shifted_h.intensity .* phase;
    end
    result = struct('intensity', h_sum, 'x', plane1.x, 'y', plane1.y);
end