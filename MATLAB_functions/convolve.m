function res_struct = convolve(im1, im2, normalize)
    %code to compute the colvoluion of 2 image structures
    arguments
        im1
        im2
        normalize = false
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
    assert(strcmp(field_type2, field_type1))
    assert(all(size(im1.(field_type1))==size(im2.(field_type2))))
    %collect image sizes to zero-pad
    im_sz = size(im1.(field_type1));
    ft_sz = im_sz;
    %fft (with padding)
    ft1 = fftshift(fft2(ifftshift(im1.(field_type1)), ft_sz(1), ft_sz(2)));
    ft2 = fftshift(fft2(ifftshift(im2.(field_type2)), ft_sz(1), ft_sz(2)));
    %multiply in fourier space
    Ft_prod = ft1 .* ft2;
    %take ift for result (and unpad)
    res = fftshift(ifft2(ifftshift(Ft_prod), im_sz(1), im_sz(2)));
    if normalize
        norm = sum(res, 'all');
        res = res ./ norm;
    end
    res_struct = struct(field_type1, res, 'x', im1.x, 'y', im1.y);
end