%% Parameters
step = 3; % Step size used during reconstruction
display_step = 9; % Display every sixth image to reduce size
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
save('reconstructed_volume_9.mat', 'volume');

%% Create the figure and the slider
hFig = figure('Name', '3D Volume Visualization with Slider');
hAx = axes('Parent', hFig);

% Initial display of the middle slice
slice_index = round(num_display_slices / 2);
hImg = imagesc(volume(:, :, slice_index), 'Parent', hAx);
colormap('gray');
colorbar;
title(hAx, sprintf('Slice %d', slice_index));
axis image;

% Slider control
hSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', num_display_slices, ...
    'Value', slice_index, 'SliderStep', [1/(num_display_slices-1), 10/(num_display_slices-1)], ...
    'Units', 'normalized', 'Position', [0.1 0.01 0.8 0.05]);

% Add a listener to the slider to update the displayed slice
addlistener(hSlider, 'Value', 'PostSet', @(src, event) updateSlice(round(get(hSlider, 'Value')), volume, hImg, hAx));

function updateSlice(slice_index, volume, hImg, hAx)
    % Update the image data and title based on the slider value
    set(hImg, 'CData', volume(:, :, slice_index));
    title(hAx, sprintf('Slice %d', slice_index));
end
