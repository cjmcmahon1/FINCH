function plot_im(image_struct, label, type, cax)
    %{
    Quick function to plot a square B&W intensity plot. Supports plotting
    of structures with 'field' or 'intensity' argument. image_struct must
    also have an 'x' and 'y' argument for labeling the axes. Type can be
    either 'intensity', 'real', or 'imag'; to get different components of
    the image plane as needed.
    %}
    arguments
        image_struct %struct with an intensity plot
        label %title for the plot
        type = 'intensity' 
        cax = 'auto' %limits for the colorbar (form is [min max])
    end
    if isfield(image_struct, 'intensity')
        field_type = 'intensity';
    elseif isfield(image_struct, 'field')
        field_type = 'field';
    else
        fprintf("Struct did not have an 'intensity' or 'field' field");
    end
    if strcmp(type,'intensity') && strcmp(field_type,'intensity')
        plot_param = abs(image_struct.(field_type));
    elseif strcmp(type,'intensity') && strcmp(field_type,'field')
        plot_param = abs(image_struct.(field_type)).^2;
    elseif strcmp(type,'real')
        plot_param = real(image_struct.(field_type));
    elseif strcmp(type,'imag')
        plot_param = imag(image_struct.(field_type));
    else
        fprintf("Did not recognize plot type. Must be 'intensity', " + ...
            "'real, or 'imag'.");
    end
    imagesc(image_struct.x*1e3, image_struct.y*1e3, ...
        plot_param);
    title(label);
    axis('square');
    colormap('gray');
    caxis(cax);
    xlabel("x (um)");
    ylabel("y (um)");
    colorbar();
end
