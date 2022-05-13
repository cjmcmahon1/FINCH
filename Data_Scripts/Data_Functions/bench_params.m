function PARAMS = bench_params(num_x_pixels, num_y_pixels, NA)
    arguments
        num_x_pixels
        num_y_pixels
        NA = 0.1
    end
    %get the bench parameters for the camera currently in use. Specify
    %number of pixels in each dimension and it will rescale it to the
    %appropriate units.
    %Measurements from the camera (units mm)
    dx = 3.45e-3; %mm
    dy = dx;
%     num_x_pixels = 1080;
%     num_y_pixels = 1440;
    L_x = dx * num_x_pixels;
    L_y = dy * num_y_pixels;
    X = -L_x/2:dx:L_x/2-dx;
    Y = -L_y/2:dy:L_y/2-dy;

    PARAMS = struct;
    PARAMS.Lx = L_x;      %x side length of input image
    PARAMS.Ly = L_y;      %y side length of input image
    PARAMS.lambda = 490e-6; %wavelength
    PARAMS.Mx = num_x_pixels;        %x samples
    PARAMS.My = num_y_pixels;        %y samples
    PARAMS.NA = NA;        %numerical aperture
    PARAMS.x = X;
    PARAMS.y = Y;
end