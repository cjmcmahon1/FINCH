function F = hologram3D_to_frames(hologram_struct, label)
    %Convert a 3D hologram into a Frames object, which can be used to
    %generate movies/scans of the resulting fresnel propagated image over
    %various z distances. (need to add z label on each frame for clarity)
    z_vals = hologram_struct.z;
    x_vals = hologram_struct.x;
    hol_max = max(abs(hologram_struct.intensity(:)));
    F(length(z_vals)) = struct('cdata',[],'colormap',[]);
    for i = 1:length(z_vals)
        X = abs(squeeze(hologram_struct.intensity(:,:,i)));
        imagesc(hologram_struct.x, hologram_struct.y, X./hol_max);
        colormap('gray');
        caxis([0 1]);
        axis('square');
        xlabel('x (mm)');
        ylabel('y (mm)');
        z_label = sprintf('z=%.2f', z_vals(i));
        title(label);
        text(x_vals(12), x_vals(12), z_label, 'BackgroundColor', 'white');
        colorbar();
        drawnow
        F(i) = getframe(gcf);
    end
end