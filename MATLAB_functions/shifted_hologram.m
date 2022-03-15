function result = shifted_hologram(plane, theta)
    %{
    Based on Brooker (2021) equation 2. Generate a hologram with a phase
    shift e^{i \theta}
    %}
%     arguments
%         plane %interference plane we get from propagate()
%         theta %artificial phase shift of the interference
%     end
    intensity = abs(plane.field).^2 .* exp(1i * theta);
    result = struct('intensity', intensity, 'x', plane.x, 'y', plane.y);
end