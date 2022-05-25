function result = shifted_hologram(plane1, plane2, theta, lambda, normalize)
    %{
    Based on Brooker (2021) equation 2. Generate a hologram with a phase
    shift e^{i \theta}
    %}
    arguments
        plane1 %interference plane we get from propagate()
        plane2
        theta %artificial phase shift of the interference
        lambda = 0;
        normalize = false;
    end
    interference  = (plane1.field .* exp(1i * theta)) + plane2.field;
    intensity = abs(interference).^2;
    if normalize %normalize the hologram to have total intensity 1
        norm = sum(intensity, 'all');
        intensity = intensity ./ norm;
    end
    if lambda > 0 %add poisson noise if a noise level is specified
        noise = poissrnd(intensity) .* lambda;
        intensity = intensity + noise;
    end
    result = struct('intensity', intensity, 'x', plane1.x, 'y', plane1.y);
end