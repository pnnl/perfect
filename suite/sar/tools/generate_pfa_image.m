% [img] = generate_pfa_image(kernel2_output_filename, image_size)
%
% Generates and displays a complex image by applying a 2D IFFT to
% the output of SAR kernel 2 from the TAV PERFECT benchmark suite.
% The arguments are as follows:
%
%    kernel2_output_filename : name of the file containing the
%        output from SAR kernel 2.
%
%    image_size : size of the final resampled data from kernel 2,
%        which will also be the size of the generated image.  This
%        should be set consistently with [PFA_NOUT_RANGE, PFA_NOUT_AZIMUTH]
%        from common/sar_params.h.
%
% This function returns the complex image.
function [img] = generate_pfa_image(kernel2_output_filename, image_size)
    fp = fopen(kernel2_output_filename, 'rb');
    img_fourier = fread(fp, inf, 'single');
    fclose(fp);

    if numel(img_fourier) ~= 2*prod(image_size)
        error('The input image size is inconsistent with the kernel output file size.');
    end

    img_fourier = complex(img_fourier(1:2:end),img_fourier(2:2:end));
    img_fourier = reshape(img_fourier, image_size);

    img = fftshift(ifft2(ifftshift(img_fourier)));

    img = flipud(fliplr(img.'));

    figure; imagesc(20*log10(abs(img./max(abs(img(:))))));
    axis xy;
    colorbar;
    caxis([-80 0]);
    title('Pixel magnitude in normalized decibel units')

end
