function result = shifted_hologram(plane, theta, bench_params, rh)
    %{
    Based on Brooker (2021) equation 2. Generate a real-valued hologram
    with a given phase shift theta from an input interference image     
    intensity. We assume the input image is of the form 
    ~exp[i/z(x.^2 + y.^2)].
    %}
%     arguments
%         plane %interference plane we get from propagate()
%         theta %artificial phase shift of the interference
%         bench_params
%         rh = 250e-3 %maximum radius of the hologram
%     end
    if isfield(plane, 'intensity')
        field_type = 'intensity';
    elseif isfield(plane, 'field')
        field_type = 'field';
    else
        fprintf("Struct did not have an 'intensity' or 'field' field");
    end
    P = pupil_func(rh, bench_params);
    h1 = plane.(field_type) .* exp(1i * theta);
    h2 = conj(plane.(field_type)) .* exp(-1i * theta);
%     intensity = P .* (2 + h1 + h2); % MS - why add 2? also, what is P doing?
    intensity = abs(h1 + h2).^2; % is this the correct expression?
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