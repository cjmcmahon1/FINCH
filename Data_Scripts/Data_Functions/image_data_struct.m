function hol = image_data_struct(image, angle)
    %create an image structure from an input intensity array. Specify the
    %relative angle for hologram reconstruction
    im_shape = size(image);
    crop_params = bench_params(im_shape(1), im_shape(2));
    hol = struct('intensity', image, 'x', crop_params.x, ...
        'y', crop_params.y, 'angle', angle);
end