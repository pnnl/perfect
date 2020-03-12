% function [img] = read_wami_image(name)
%
% Reads WAMI image files that include a header with the
% [width height channel depth] as the first four uint16 fields
% in the file.  The input and output images for kernel 1 (debayer)
% follow this format.
function [img] = read_wami_image(name)

    fp = fopen(name, 'rb');
    width = fread(fp, 1, 'uint16');
    height = fread(fp, 1, 'uint16');
    chan = fread(fp, 1, 'uint16');
    bytes_per_pixel = fread(fp, 1, 'uint16');

    if bytes_per_pixel ~= 2, error('Unhandled bit depth'); end

    data = fread(fp, inf, 'uint16');
    fclose(fp);

    img = reshape(data, [chan width height]);
    img = permute(img, [3 2 1]);
    figure;
    imagesc(img ./ max(img(:)));
    if chan == 1, colormap gray; end

end
