clc;
clear;

%% Load the saved sinogram
load('sinogram.mat', 'sinogram');

% Number of angles (images)
num_angles = 400; % 400 images
theta = linspace(0, 360, num_angles); % 400 angles in degrees
theta_rad = theta * pi / 180;  % Convert angles to radians

% Normalize the sinogram for better visualization
sinogram_normalized = mat2gray(sinogram);

% Visualize the sinogram
figure;
imshow(sinogram_normalized, []), title('Normalized Sinogram at Height 1120');
colormap(gca, 'gray'); % Set the colormap to gray for better visualization
colorbar; % Add a colorbar for reference