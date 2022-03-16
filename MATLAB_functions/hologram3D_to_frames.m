function F = hologram3D_to_frames(hologram_struct)
    %Convert a 3D hologram into a Frames object, which can be used to
    %generate movies/scans of the resulting fresnel propagated image over
    %various z distances. (need to add z label on each frame for clarity)
    z_vals = hologram_struct.z;
    hol_max = max(abs(hologram_struct.intensity(:)));
    F(length(z_vals)) = struct('cdata',[],'colormap',[]);
    for i = 1:length(z_vals)
        X = abs(squeeze(hologram_struct.intensity(:,:,i)));
        imagesc(hologram_struct.x, hologram_struct.y, X./hol_max);
        colormap('gray');
        caxis([0 1]);
        axis('square');
        xlabel("x (mm)");
        ylabel("y (mm)");
        colorbar();
        drawnow
        F(i) = getframe;
    end
end