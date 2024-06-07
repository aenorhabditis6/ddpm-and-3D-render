function cone_beam_metal_sino_filter(filter)
    % Parameters
    num_angles = 1081; 
    height = 512; % Height of the images
    middle_slice_index = 256; % Middle slice index
    D1 = 55.8595; % Source to object distance
    D2 = 1115.956 - 55.8595; % Object to detector distance

    % Read the projection images
    projections = zeros(height, height, num_angles);
    for i = 1:num_angles
        filename = sprintf('C:/Users/tobeo/Desktop/metal_512/Cylinder_0.1mmAl_150kV_90uA_2s_24db_32_5_1080_4_RAC_0mm_%04d.tif', i);
        projections(:, :, i) = imread(filename);
    end

    % Extract the cross-section and build the sinogram
    sino = zeros(height, num_angles); % Initialize sinogram with dimensions swapped
    for i = 1:num_angles
        % Extract the cross-section
        sino(:, i) = projections(middle_slice_index, :, i);
    end

    % Save the original sinogram
    sinogram_filename = './sinogram_middle_slice_256_original.tif';
    imwrite(mat2gray(sino), sinogram_filename);
    fprintf('Original Sinogram saved at %s.\n', sinogram_filename);

    % Frequency scaling
    d = D1 + D2;
    w = (0:size(sino,1)-1)' * pi / size(sino,1);
    filt = abs(w);

    % Apply selected filter
    switch lower(filter)
        case 'ram-lak'
            % Do nothing
        case 'shepp-logan'
            % be careful not to divide by 0:
            filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
        case 'cosine'
            filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
        case 'hamming'  
            filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
        case 'hann'
            filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
        otherwise
            error('Invalid filter selected.');
    end

    % Apply the filter to the sinogram
    fft_sino = fft(sino);
    filt = repmat(filt, 1, size(fft_sino, 2));
    filtered_sino = real(ifft(fft_sino .* filt));

    % Save the filtered sinogram
    filtered_sinogram_filename = sprintf('./sinogram_middle_slice_256_%s.tif', filter);
    imwrite(mat2gray(filtered_sino), filtered_sinogram_filename);
    fprintf('Filtered Sinogram (%s) saved at %s.\n', filter, filtered_sinogram_filename);

    % (Reconstruction steps would go here using 'filtered_sino')
end
