function result = mirror_hologram(plane1, plane2, num_angles)
    %{
    Identical to complex_hologram.m but with a flipped sign in the phase.
    This generates the twin image for comparison.
    %}
    arguments
        plane1 %interference plane we get from propagate()
        plane2 %second interference plane we get from propagate ()
        num_angles %>=3 subdivisions of 2*pi to apply differing phases
    end
    inc = 2*pi / num_angles; %we want (num_angles) evenly separated phases
    h_sum = zeros('like', plane1.field);
    for i = 1:num_angles
        prev_angle = mod(inc*(i-1), 2*pi);
        next_angle = mod(inc*(i+1), 2*pi);
        shifted_h = shifted_hologram(plane1, plane2, inc*i);
        phase = exp(-1i * prev_angle) - exp(-1i * next_angle); %sign flip
        h_sum = h_sum + shifted_h.intensity .* phase;
    end
    result = struct('intensity', h_sum, 'x', plane1.x, 'y', plane1.y);
end