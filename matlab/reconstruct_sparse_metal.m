clc;
clear;
num_angles = 1081; 
height = 2000; 
step = 2; 
num_images = ceil(num_angles / step); % Number of images to process
projections = zeros(2000, 2000, num_images); 

for i = 1:step:num_angles
    index = ceil(i / step); 
    filename = sprintf('./AM/Cylinder_0.1mmAl_150kV_90uA_2s_24db_32_5_1080_4_RAC_0mm_%04d.tif', i);
    projections(:, :, index) = imread(filename);
    fprintf('Processing file: %s\n', filename);
end

slice_counter = 1;
for slice_index = 1:step:height
    fprintf('Processing slice %d/%d...\n', slice_index, height);
    
    sino = zeros(2000, num_images); 
    for i = 1:num_images
        sino(:, i) = projections(slice_index, :, i);
    end

    %% Perform Radon Transform
    R = sino;
    width = 2^nextpow2(size(R, 1));

    %% Perform FFT and filtering on the projections
    proj_fft = fft(R, width);

    % Ramp filter
    filter = 2 * [0:(width / 2 - 1), width / 2:-1:1]' / width;

    % Filtered projections
    proj_filtered = zeros(width, num_images);
    for i = 1:num_images
        proj_filtered(:, i) = proj_fft(:, i) .* filter;
    end

    %% IFFT
    proj_ifft = real(ifft(proj_filtered));

    %% Backprojection
    M = 2000; 
    fbp = zeros(M);
    theta_rad = linspace(0, 360, num_images) * pi / 180; % Convert angles to radians
    xp_offset = ceil(size(proj_ifft, 1) / 2);

    for i = 1:num_images
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
    fbp = (fbp * pi) / num_images;

    fbp_normalized = mat2gray(fbp);
   
    reconstructed_image_filename = sprintf('./reconstructed_metal/reconstructed_image_%d_ramp.png', slice_index);
    imwrite(fbp_normalized, reconstructed_image_filename);
end

fprintf('Reconstruction complete.\n');
