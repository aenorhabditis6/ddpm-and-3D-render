function reconstruct_full_volume()
    % Parameters
    num_angles = 400; % 400 images
    height = 1952; % Height of the images
    
    % Read the projection images
    projections = zeros(1952, 1952, num_angles);
    for i = 1:num_angles
        filename = sprintf('./TIF/g_%03d.tif', i);
        projections(:, :, i) = imread(filename);
    end
    
    % Initialize 3D volume
    volume = zeros(1952, 1952, height);
    
    % Create the directory for saving reconstructed images if it doesn't exist
    if ~exist('./reconstructed', 'dir')
        mkdir('./reconstructed');
    end
    
    % Loop over each slice index
    for slice_index = 1:height
        fprintf('Processing slice %d/%d...\n', slice_index, height);
        
        %% Extract the cross-section and build the sinogram
        sino = zeros(1952, num_angles); % Initialize the sinogram
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

        %% Backprojection
        M = 1952; % Number of pixels in the reconstructed image
        fbp = zeros(M);
        theta_rad = linspace(0, 360, num_angles) * pi / 180;  % Convert angles to radians
        xp_offset = ceil(size(proj_ifft, 1) / 2);

        for i = 1:num_angles
            rad = theta_rad(i); % Angle in radians
            for x = 1:M
                for y = 1:M
                    t = round((x - M/2) * cos(rad + pi/2) - (y - M/2) * sin(rad + pi/2) + xp_offset);
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
        reconstructed_image_filename = sprintf('./reconstructed/reconstructed_image_%d_ramp.png', slice_index);
        imwrite(fbp_normalized, reconstructed_image_filename);

        % Add to the 3D volume
        volume(:, :, slice_index) = fbp_normalized;
    end
    
    % Save the 3D volume
    save('reconstructed_volume.mat', 'volume');
    
    fprintf('Reconstruction complete. 3D volume saved as reconstructed_volume.mat\n');
end