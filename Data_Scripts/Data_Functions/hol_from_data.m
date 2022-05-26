function c_hologram = hol_from_data(image_struct_list)
    % create a complex valued hologram from an input array of
    % image_data_struct() structures. The structures must have an 'angle'
    % field to specify the relative angles between them.
    num_h = length(image_struct_list);
    int_shape = size(image_struct_list(1).intensity);
    angles = zeros(1, num_h);
    intensities = zeros(num_h, int_shape(1), int_shape(2));
    for i = 1:num_h
        angles(i) = image_struct_list(i).angle;
        intensities(i,:,:) = [image_struct_list(i).intensity];
    end
    h_sum = zeros(int_shape);
    for i = 1:num_h
        prev_angle_idx = mod(i - 2 + num_h, num_h)+1;
        prev_angle = angles(prev_angle_idx);
        next_angle_idx = mod(i, num_h)+1;
        next_angle = angles(next_angle_idx);
        phase = exp(1i * prev_angle) - exp(1i * next_angle);
        h_sum = h_sum + squeeze(intensities(i,:,:)) .* phase;
    end
    c_hologram = struct('intensity', h_sum, 'x', ...
        image_struct_list(1).x, 'y', image_struct_list(1).y);
end