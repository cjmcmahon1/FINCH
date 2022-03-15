function hol_struct = gen_hol_im(hol, z, bench_params)
    % generate a structure by fresnel propagating an input hologram.
    % Parameters are:
    % hol: complex valued hologram created through hol_from_data(), or
    % complex_hologram()
    % z: distance to z propagate hologram
    % bench_params: bench parameters
    hol_plane = fresnel_prop(hol.intensity, z, bench_params);
    hol_struct = struct('intensity', hol_plane, 'x', hol.x, 'y', hol.y);
end