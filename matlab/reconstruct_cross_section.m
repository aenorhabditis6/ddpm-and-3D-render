function reconstruct_cross_section(slice_index)
    %% Parameters
    num_angles = 400; % 400 images

    % Read the projection images
    projections = zeros(1952, 1952, num_angles);
    for i = 1:num_angles
        filename = sprintf('./TIF/g_%03d.tif', i);
        projections(:, :, i) = imread(filename);
    end

    %% Extract the cross-section and build the sinogram
    sino = zeros(1952, num_angles); % Initialize the sinogram
    for i = 1:num_angles
        % Extract the cross-section
        sino(:, i) = projections(slice_index, :, i);
    end

    % Normalize the sinogram for better visualization
    sinogram_normalized = mat2gray(sino);

    % Visualize the sinogram
    figure;
    imshow(sinogram_normalized, []), title(sprintf('Normalized Sinogram at Height %d', slice_index));
    colormap(gca, 'gray');
    colorbar; % See normalized value

    %% Save the sinogram for future use
    sinogram_filename = sprintf('sinogram_%d.mat', slice_index);
    save(sinogram_filename, 'sino');

    %% Load the saved sinogram
    load(sinogram_filename, 'sino');

    theta = linspace(0, 360, num_angles);
    theta_rad = theta * pi / 180;  % Convert angles to radians

    % Parameters for reconstruction
    M = 1952; % Number of pixels in the reconstructed image

    % Perform Radon Transform
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

    figure;
    subplot(1, 2, 1), imshow(log(abs(proj_fft) + 1), []), title('FFT of Projections');
    subplot(1, 2, 2), imshow(log(abs(proj_filtered) + 1), []), title('Filtered FFT Projections');

    %% Perform inverse FFT
    proj_ifft = real(ifft(proj_filtered));

    figure;
    imshow(proj_ifft, []), title('Inverse FFT of Filtered Projections');

    %% Backprojection
    fbp = zeros(M);
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

    % Display the reconstructed image
    figure;
    imshow(fbp, []), title(sprintf('Reconstructed Image at Height %d', slice_index));

    %% Save the reconstructed image
    reconstructed_image_filename = sprintf('reconstructed_image_%d.png', slice_index);
    imwrite(mat2gray(fbp), reconstructed_image_filename);
    
    % excel_filename = sprintf('rebuild_info_%d.xlsx', slice_index);
    % xlswrite(excel_filename, fbp, 'Sheet1');
end