function result = shifted_hologram(plane1, plane2, theta, lambda)
    %{
    Based on Brooker (2021) equation 2. Generate a hologram with a phase
    shift e^{i \theta}
    %}
    arguments
        plane1 %interference plane we get from propagate()
        plane2
        theta %artificial phase shift of the interference
        lambda = 0;
    end
    interference  = (plane1.field .* exp(1i * theta)) + plane2.field;
    intensity = abs(interference).^2;
    if lambda > 0 %add poisson noise if a noise level is specified
        noise = poissrnd(lambda, size(intensity, 1), size(intensity, 2));
        intensity = intensity + noise;
    end
    result = struct('intensity', intensity, 'x', plane1.x, 'y', plane1.y);
end