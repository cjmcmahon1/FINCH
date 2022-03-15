function result = shifted_hologram(plane1, plane2, theta)
    %{
    Based on Brooker (2021) equation 2. Generate a hologram with a phase
    shift e^{i \theta}
    %}
%     arguments
%         plane %interference plane we get from propagate()
%         theta %artificial phase shift of the interference
%     end
    interference  = (plane1.field .* exp(1i * theta)) + plane2.field;
    intensity = abs(interference).^2;
    result = struct('intensity', intensity, 'x', plane1.x, 'y', plane1.y);
end