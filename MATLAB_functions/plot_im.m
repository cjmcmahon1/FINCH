function plot_im(image_struct, label)
    %quick function to plot a square B&W intensity plot
    imagesc(image_struct.x*1e3, image_struct.y*1e3, ...
        abs(image_struct.field).^2);
    title(label);
    axis('square');
    colormap('gray');
    xlabel("x (um)");
    ylabel("y (um)");
    colorbar();
end
