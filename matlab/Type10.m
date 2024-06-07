% filter_econstruct_cross_section(1320, 'hann');
% Run the function
reconstruct_sparse_volume();

%%
load('reconstructed_volume.mat', 'volume');
volshow(volume);

%%
clc;
clear;

%% Parameters
step = 3; % Step size used during reconstruction
display_step = 12; % Display every sixth image to reduce size
height = 1952; % Height of the original images
num_slices = ceil(height / step); % Number of slices reconstructed
num_display_slices = ceil(num_slices / display_step); % Number of slices to display

% Initialize 3D volume with reduced size
volume = zeros(1952, 1952, num_display_slices);

% Read each reconstructed image and add to the 3D volume
slice_counter = 1;
for slice_index = 1:display_step*step:height
    fprintf('Reading slice %d...\n', slice_index);
    
    % Construct the filename
    reconstructed_image_filename = sprintf('./reconstructed/reconstructed_image_%d_ramp.png', slice_index);
    
    % Read the image
    img = imread(reconstructed_image_filename);
    
    % Convert to double precision and normalize if needed
    img = im2double(img);
    
    % Add to the 3D volume
    volume(:, :, slice_counter) = img;
    slice_counter = slice_counter + 1;
end

% Save the 3D volume
save('reconstructed_volume_12.mat', 'volume');

%% Visualize the 3D volume
if exist('volshow', 'file')
    volshow(volume);
else
    % Fall back to isosurface visualization
    figure;
    p = patch(isosurface(volume, 0.5));
    isonormals(volume, p);
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1]);
    view(3); 
    axis tight;
    camlight; 
    lighting gouraud;
    title('Isosurface Visualization of the Reconstructed 3D Volume');
end


%% 

%% Visualize the 3D volume using orthogonal slices
figure;
colormap('gray');

% Create slices along the x, y, and z axes
slice(double(volume), size(volume, 2)/2, size(volume, 1)/2, size(volume, 3)/2);
shading interp;
colorbar;
title('Orthogonal Slice Visualization of the Reconstructed 3D Volume');

% Adjust view settings for better visualization
daspect([1 1 1]);
view(3);
axis tight;
camlight;
lighting gouraud;

%%

im = imread('./AM/Cylinder_0.1mmAl_150kV_90uA_2s_24db_32_5_1080_4_RAC_0mm_0001.tif');
