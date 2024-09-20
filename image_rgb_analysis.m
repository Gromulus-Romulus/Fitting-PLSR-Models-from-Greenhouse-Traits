% Objective: Calculate
% green, % yellow, % brown
% of leaf scan (whole images).
%
% Date: 2024.05.25
% Author: Nathan
% Written with ChatGPT v4
%
% Source: 
%   https://chatgpt.com/share/6799c940-667d-44a5-bfe0-123880cacb76
%

% Get a list of all image files in the scans folder
folder = "./scans/";
imageFiles = dir(fullfile(folder, '*.jpg')); % Adjust the extension if needed

% Main Loop
for k = 1:length(imageFiles)
    % Read the image
    img = imread(fullfile(folder, imageFiles(k).name));
    
    % Convert the image to the HSV color space for easier color segmentation
    hsvImg = rgb2hsv(img);
    
    % Define thresholds for green and yellow colors in HSV
    greenThresholdLow = [0.25, 0.20, 0.20];  % Adjust as needed
    greenThresholdHigh = [0.50, 1.00, 1.00]; % Adjust as needed
    yellowThresholdLow = [0.10, 0.20, 0.20]; % Adjust as needed
    yellowThresholdHigh = [0.20, 1.00, 1.00]; % Adjust as needed
    
    % Create binary masks for green and yellow regions
    greenMask = (hsvImg(:,:,1) >= greenThresholdLow(1)) & (hsvImg(:,:,1) <= greenThresholdHigh(1)) & ...
                (hsvImg(:,:,2) >= greenThresholdLow(2)) & (hsvImg(:,:,2) <= greenThresholdHigh(2)) & ...
                (hsvImg(:,:,3) >= greenThresholdLow(3)) & (hsvImg(:,:,3) <= greenThresholdHigh(3));
            
    yellowMask = (hsvImg(:,:,1) >= yellowThresholdLow(1)) & (hsvImg(:,:,1) <= yellowThresholdHigh(1)) & ...
                 (hsvImg(:,:,2) >= yellowThresholdLow(2)) & (hsvImg(:,:,2) <= yellowThresholdHigh(2)) & ...
                 (hsvImg(:,:,3) >= yellowThresholdLow(3)) & (hsvImg(:,:,3) <= yellowThresholdHigh(3));
    
    % Calculate the number of pixels for each color
    totalPixels = numel(hsvImg(:,:,1));
    greenPixels = sum(greenMask(:));
    yellowPixels = sum(yellowMask(:));
    
    % Calculate the percentage of each color
    percentGreen = (greenPixels / totalPixels) * 100;
    percentYellow = (yellowPixels / totalPixels) * 100;
    
    % Display the results
    fprintf('Image %s:\n', imageFiles(k).name);
    fprintf('Percentage of Green: %.2f%%\n', percentGreen);
    fprintf('Percentage of Yellow: %.2f%%\n', percentYellow);

end;
