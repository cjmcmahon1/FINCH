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
    %intensity = P .* (2 + h1 + h2); % MS - why add 2? also, what is P doing?
    result = struct('intensity', intensity, 'x', plane.x, 'y', plane.y);
end

function plane = pupil_func(radius, bench_params)
    %pupil function in real space
%     arguments
%         radius %mm
%         bench_params
%     end
    dx = bench_params.Lx/bench_params.Mx;
    x = -bench_params.Lx/2:dx:bench_params.Lx/2-dx;
    dy = bench_params.Ly/bench_params.My;
    y = -bench_params.Ly/2:dy:bench_params.Ly/2-dy;
    [X,Y] = meshgrid(x,y);
    plane = (X.^2 + Y.^2) <= radius^2;
end