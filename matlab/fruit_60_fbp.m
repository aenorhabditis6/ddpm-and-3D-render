function fruit_60_fbp()
    % Parameters
    num_angles = 133; % 400 / 3 images
    height = 512; % Height of the images
    step = 5; % Step
    D1 = 1573.724; % Source to object distance
    D2 = 399.899; % Object to detector distance

    % Read the projection images
    projections = zeros(height, height, num_angles);
    for i = 1:num_angles
        filename = sprintf('C:/Users/tobeo/Desktop/TIF_512/g_%03d.tif', i);
        projections(:, :, i) = imread(filename);
    end

    % Calculate the number of slices to be reconstructed
    num_slices = ceil(height / step);

    % Create the directory for saving reconstructed images if it doesn't exist
    if ~exist('./partially_reconstruct_fruit_fbp', 'dir')
        mkdir('./partially_reconstruct_fruit_fbp');
    end

    % Precompute the angles in radians and their sines and cosines
    theta_rad = linspace(0, 120, num_angles) * pi / 180;
    cos_theta = cos(theta_rad);
    sin_theta = sin(theta_rad);

    for slice_index = 1:step:height
        fprintf('Processing slice %d/%d...\n', slice_index, height);
        
        %% Extract the cross-section and build the sinogram
        sino = zeros(512, num_angles); % Initialize the sinogram
        for i = 1:num_angles
            % Extract the cross-section
            sino(:, i) = projections(slice_index, :, i);
        end

        %% Perform Radon Transform
        R = sino;

        % Set the width for the Fast Fourier Transform (FFT)
        width = 2^nextpow2(size(R, 1));

        %% Perform FFT and filtering on the projections
        proj_fft = fft(R, width);

        % Ramp filter
        filter = 2 * [0:(width / 2 - 1), width / 2:-1:1]' / width;

        % Filtered projections
        proj_filtered = zeros(width, num_angles);
        for i = 1:num_angles
            proj_filtered(:, i) = proj_fft(:, i) .* filter;
        end

        %% IFFT
        proj_ifft = real(ifft(proj_filtered));

        %% Backprojection for Cone-Beam Geometry
        M = 512; % Number of pixels in the reconstructed image
        fbp = zeros(M);
        xp_offset = ceil(size(proj_ifft, 1) / 2);
        
        for i = 1:num_angles
            rad = theta_rad(i); % Angle in radians
            cos_t = cos_theta(i);
            sin_t = sin_theta(i);
            for x = 1:M
                for y = 1:M
                    % Calculate the t value considering cone-beam geometry
                    z = slice_index - height / 2; % Adjust for slice index
                    X = (x - M / 2) * D1 / (D1 + z);
                    Y = (y - M / 2) * D1 / (D1 + z);
                    t = round((X * cos_t + Y * sin_t) + xp_offset);
                    if t > 0 && t <= size(proj_ifft, 1)
                        fbp(x, y) = fbp(x, y) + proj_ifft(t, i);
                    end
                end
            end
        end
        fbp = (fbp * pi) / num_angles;

        % Normalize the reconstructed image
        fbp_normalized = mat2gray(fbp);
        
        % Save the reconstructed image
        reconstructed_image_filename = sprintf('./partially_reconstruct_fruit_fbp/reconstructed_512_fruit_%d_120.png', slice_index);
        imwrite(fbp_normalized, reconstructed_image_filename);

        % Add to the 3D volume
        % volume(:, :, slice_counter) = fbp_normalized;
        % slice_counter = slice_counter + 1;
    end
    
    % Save the 3D volume
    % save('reconstructed_volume.mat', 'volume');
    % 
    fprintf('Reconstruction complete.');
end
