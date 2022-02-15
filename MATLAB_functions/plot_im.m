function plot_im(image_struct, label)
    %quick function to plot a square B&W intensity plot
    if isfield(image_struct, 'intensity')
        plot_param = abs(image_struct.intensity);
    elseif isfield(image_struct, 'field')
        plot_param = abs(image_struct.field).^2;
    end
    imagesc(image_struct.x*1e3, image_struct.y*1e3, ...
        plot_param);
    title(label);
    axis('square');
    colormap('gray');
    xlabel("x (um)");
    ylabel("y (um)");
    colorbar();
end
