function mat = open_im(filename)
    %convert a PNG into a double (0-1) array (grayscale).
    im = imread(filename);
    mat = im2double(rgb2gray(im));
end