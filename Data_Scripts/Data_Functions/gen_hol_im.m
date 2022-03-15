function hol_struct = gen_hol_im(hol, z, bench_params)
    % generate a hologram image propagated a given distance z
    hol_plane = fresnel_prop(hol.intensity, z, bench_params);
    hol_struct = struct('intensity', hol_plane, 'x', hol.x, 'y', hol.y);
end