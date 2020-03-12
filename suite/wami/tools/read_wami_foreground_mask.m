% function [img] = read_wami_foreground_mask(name, paramset)
%
% Reads the foreground mask generated as output for WAMI kernel 3
% (change detection / GMM) from the PERFECT suite.  name is
% the filename to read and paramset should be 'small', 'medium',
% or 'large', depending on the parameter set used.  The output
% img volume will have dimension [height width nframes].
function [img] = read_wami_foreground_mask(name, paramset)

    if nargin ~= 2
        fprintf(1, 'Invalid input parameters; help follows\n');
        eval(sprintf('help %s', mfilename));
        return;
    end

    fp = fopen(name, 'rb');
    data = fread(fp, inf, 'uint8');
    fclose(fp);

    if strcmpi(paramset, 'small')
        width = 512;
        height = 512;
    elseif strcmpi(paramset, 'medium')
        width = 1024;
        height = 1024;
    elseif strcmpi(paramset, 'large')
        width = 2048;
        height = 2048;
    else
        fprintf(1, 'Error: unsupported paramset; should be small, medium, or large.\n');
        return;
    end

    nframes = numel(data) / (width * height);
    img = reshape(data, [width height nframes]);
    for i = 1:nframes, img(:,:,i) = img(:,:,i)'; end
    figure;
    imagesc(img(:,:,1));
    colormap gray;

end
