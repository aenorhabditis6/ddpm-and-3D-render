function cone_beam_metal_2()
    % Parameters
    num_angles = 1081; 
    height = 512; % Height of the images
    step = 5; % Step
    D1 = 55.8595; % Source to object distance
    D2 = 1115.956 - 55.8595; % Object to detector distance
%     SrcToObject=55.8595
% SrcToDetector=1115.956

    % Read the projection images
    projections = zeros(height, height, num_angles);
    for i = 1:num_angles
        filename = sprintf('C:/Users/tobeo/Desktop/metal_512/Cylinder_0.1mmAl_150kV_90uA_2s_24db_32_5_1080_4_RAC_0mm_%04d.tif', i);
        projections(:, :, i) = imread(filename);
    end

    % Create the directory for saving reconstructed images if it doesn't exist
    if ~exist('./metal_512_fix_5_2', 'dir')
        mkdir('./metal_512_fix_5_2');
    end

    % Precompute the angles in radians and their sines and cosines
    theta_rad = linspace(0, 359, num_angles) * pi / 180;
    cos_theta = cos(theta_rad);
    sin_theta = sin(theta_rad);

    % Loop over each slice index with the specified step
    for slice_index = 1:step:height
        fprintf('Processing slice %d/%d...\n', slice_index, height);
        
        % Extract the cross-section and build the sinogram
        sino = zeros(512, num_angles);
        for i = 1:num_angles
            % Extract the cross-section
            sino(:, i) = projections(slice_index, :, i);
        end

        % Perform Radon Transform
        R = sino;

        % Set the width for the Fast Fourier Transform (FFT)
        width = 2^nextpow2(size(R, 1));

        % Perform FFT and filtering on the projections
        proj_fft = fft(R, width);

        % Ramp filter
        filter = 2 * [0:(width / 2 - 1), width / 2:-1:1]' / width;

        % Filtered projections
        proj_filtered = zeros(width, num_angles);
        for i = 1:num_angles
            proj_filtered(:, i) = proj_fft(:, i) .* filter;
        end

        % IFFT
        proj_ifft = real(ifft(proj_filtered));

        % Backprojection for Cone-Beam Geometry
        M = 512; % Number of pixels in the reconstructed image
        fbp = zeros(M);
        xp_offset = ceil(size(proj_ifft, 1) / 2);

        % Correct scaling factor for the coordinates
        % scaling_factor = D1 / (D1 + slice_index - height / 2);
        % SrcToDetector / (SrcToDetector - SrcToObject);
        scaling_factor = D2 / (D2 + D1 - height / 2);
        
        for i = 1:num_angles
            rad = theta_rad(i); % Angle in radians
            cos_t = cos_theta(i);
            sin_t = sin_theta(i);
            for x = 1:M
                for y = 1:M
                    % Calculate the t value considering cone-beam geometry
                    
                    X = (x - M / 2) * scaling_factor;
                    Y = (y - M / 2) * scaling_factor;
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
        reconstructed_image_filename = sprintf('./metal_512_fix_5_2/Cylinder_0.1mmAl_150kV_90uA_2s_24db_32_5_1080_4_RAC_0mm_%04d_512.tif', slice_index);
        imwrite(fbp_normalized, reconstructed_image_filename);
    end
    
    fprintf('Reconstruction complete.');
end
