% Objective: Remove noise from image backgrounds
%
% Date: 2024.09.26
% Author: Nathan

% Define input and output folders

inputFolder = "./raw/scans/";
outputFolder = "./raw/masked/";
imageFiles = dir(fullfile(inputFolder, '*.jpg')); % Adjust the extension if needed

% Create the output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
% If it does exist, overwrite it
else
    rmdir(outputFolder, 's');
    mkdir(outputFolder);
end

% Initialize an empty cell array to hold the results
results = {};

% Main Loop for each image
for k = 1:length(imageFiles)
    % Read the image
    img = imread(fullfile(inputFolder, imageFiles(k).name));

    % Extract the barcode ID from the filename (assuming it's the part before the extension)
    [~, barcodeID, ~] = fileparts(imageFiles(k).name);
    
    % Convert to grayscale to detect the white background
    grayImg = rgb2gray(img);
    
    % Threshold to detect the non-white regions (leaves)
    % Assuming white background is close to 255 in intensity
    leafMask = grayImg < 200; % Adjust threshold if needed

    % Remove small noise (non-leaf fragments) using morphological operations
    leafMask = bwareaopen(leafMask, 500); % Remove small objects with fewer than 500 pixels
    
    % Apply the mask to create an ROI for the leaves
    roiImg = img;
    roiImg(repmat(~leafMask, [1, 1, 3])) = 255; % Set background to white
    
    % Save the ROI image to the output folder
    [~, baseFileName, ~] = fileparts(imageFiles(k).name); % Get the base filename
    outputFileName = fullfile(outputFolder, strcat(baseFileName, '.png')); % Output file path
    imwrite(roiImg, outputFileName, FMT="PNG"); % Save the image as transparent png

end