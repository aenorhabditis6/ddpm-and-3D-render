clc;
clear;

%% Parameters
% Number of angles (images)
num_angles = 400; % 400 images
slice_index = 1320;  % The height at which I am interested in the cross-section

projections = zeros(1952, 1952, num_angles);

for i = 1:num_angles
    filename = sprintf('./TIF/g_%03d.tif', i);
    projections(:, :, i) = imread(filename);
end

%% Extract the cross-section and build the sinogram
sino = zeros(1952, num_angles); % Initialize the sinogram

for i = 1:num_angles
    % Extract the cross-section (1120th line)
    sino(:, i) = projections(slice_index, :, i);
end

sinogram_normalized = mat2gray(sino);

% Visualize the sinogram
figure;
imshow(sinogram_normalized, []), title('Normalized Sinogram at Height 1120');
colormap(gca, 'gray');
colorbar; % See normalized value, beforehand it has 10000-50000 value

%% Save the sinogram for future use
save('sinogram_1320.mat', 'sino');


%% anx
% im = imread('./TIF/g_001.tif');

%% Load the saved sinogram
load('sinogram_1320.mat', 'sino');

% Number of angles (images)
num_angles = 400; % 400 images
theta = linspace(0, 360, num_angles); % 400 angles in degrees
theta_rad = theta * pi / 180;  % Convert angles to radians

% Normalize the sinogram for better visualization
sinogram_normalized = mat2gray(sino);

% Visualize the sinogram
figure;
imshow(sinogram_normalized, []), title('Normalized Sinogram at Height 1120');
colormap(gca, 'gray'); % Set the colormap to gray for better visualization
colorbar; % Add a colorbar for reference


%% Parameters for reconstruction
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
subplot(1, 2, 2), imshow(log(abs(proj_filtered) + 1), []), title('Filtered Projections');


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
imshow(fbp, []), title('Reconstructed Image');

%% Save the reconstructed image as a PNG file
imwrite(mat2gray(fbp), 'reconstructed_image.png');

% Save the reconstructed image data as an Excel file
xlswrite('rebuild_info.xlsx', fbp, 'Sheet1');
