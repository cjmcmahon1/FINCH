function hol3d_struct = hologram3D(hologram, z_vals, bench_params)
    %generate a rank-3 complex hologram by fresnel propagating along
    %z_vals (assumed to be evenly spaced)
    num_z = length(z_vals);
    hol_shape = size(hologram.intensity);
    hol3d = zeros(hol_shape(1), hol_shape(2),num_z);
    for i = 1:num_z
        hol3d(:,:,i) = fresnel_prop(hologram.intensity, ...
                                    z_vals(i), bench_params);
    end
    hol3d_struct = struct('intensity', hol3d, 'x', hologram.x, 'y', ...
                  hologram.y, 'z', z_vals);
end