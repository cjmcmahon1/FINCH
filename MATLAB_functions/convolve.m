function res_struct = colvolve(im1, im2)
    %code to compute the colvoluion of 2 image structures
    arguments
        im1 %array that we want to propagate
        im2
    end
    if isfield(im1, 'intensity')
        field_type1 = 'intensity';
    elseif isfield(im1, 'field')
        field_type1 = 'field';
    end
    if isfield(im2, 'intensity')
        field_type2 = 'intensity';
    elseif isfield(im2, 'field')
        field_type2 = 'field';
    end
    %check array dimensions
    assert(all(size(im1.(field_type1))==size(im2.(field_type2))))
    %fft
    ft1 = fftshift(fft2(ifftshift(im1.(field_type1))));
    ft2 = fftshift(fft2(ifftshift(im2.(field_type2))));
    %multiply in fourier space
    Ft_prod = ft1 .* ft2;
    %take ift for result
    res = fftshift(ifft2(ifftshift(Ft_prod)));
    res_struct = struct(field_type1, res, 'x', im1.x, 'y', im1.y);
end