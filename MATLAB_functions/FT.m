function plane_struct = FT(image_struct)
    %Fourier Transform an image structure. Returns a structure with
    %properly scaled axes
    if isfield(image_struct, 'intensity')
        field_type = 'intensity';
    elseif isfield(image_struct, 'field')
        field_type = 'field';
    else
        fprintf("Struct did not have an 'intensity' or 'field' parameter");
    end
    ft = fftshift(fftn(ifftshift(image_struct.(field_type))));
    %get correct frequency axis
    % Define spatial axes
    dx = image_struct.x(2) - image_struct.x(1);
    dy = image_struct.y(2) - image_struct.y(1);
    lx = image_struct.x(end) - image_struct.x(1);
    ly = image_struct.y(end) - image_struct.y(1);
    % Define frequency axes
    fMax_x = 1/(2*dx);
    fMax_y = 1/(2*dy);
    df_x = 1/lx;
    df_y = 1/ly;
    fx = -fMax_x:df_x:fMax_x-df_x;
    fy = -fMax_y:df_y:fMax_y-df_y;
    plane_struct = struct('intensity', ft, 'fx', fx, 'fy', fy);
    if isfield(image_struct, 'z')
       %compute z frequencies if it's a 3D hologram
       dz = image_struct.z(2) - image_struct.z(1); 
       fMax_z = 1/(2*dz);
       lz = image_struct.z(end) - image_struct.z(1);
       df_z = 1/lz;
       fz = -fMax_z:df_z:fMax_z-df_z;
       plane_struct.fz = fz;
    end
end