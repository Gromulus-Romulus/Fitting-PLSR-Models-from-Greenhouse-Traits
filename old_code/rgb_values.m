% Objective: Output CSV spreadsheet with ONLY RGB values.
% Date: 2024.09.26
% Author: Nathan Malamud
% Written with ChatGPT4
% https://chatgpt.com/share/66f5c5a9-8f28-8012-a8aa-00f67faba522

% Define the folder containing images
folder = './masked/';

% Get a list of all image files in the folder
imageFiles = dir(fullfile(folder, '*.jpg')); % Adjust the extension if needed

% Initialize an empty cell array to hold the results
results = {};

% Loop through each image file
for k = 1:length(imageFiles)
    % Read the image
    img = imread(fullfile(folder, imageFiles(k).name));
    
    % Extract the barcode ID from the filename (assuming it's the part before the extension)
    [~, barcodeID, ~] = fileparts(imageFiles(k).name);

    % Step 1: Convert to grayscale (8-bit)
    grayImage = rgb2gray(img);

    % Step 2: Apply threshold (0 to 214)
    thresholdValue = 214;
    binaryMask = grayImage <= thresholdValue;

    % Step 3: Remove small fragments and create mask for ROI (connected components)
    % (Optional: Adjust region size threshold as needed)
    cc = bwconncomp(binaryMask);
    regionStats = regionprops(cc, 'Area');
    largeRegionsIdx = find([regionStats.Area] > 50);  % Filter small regions (adjust size)

    % Step 4: Create final binary mask (with larger regions)
    finalMask = ismember(labelmatrix(cc), largeRegionsIdx);

    % Step 5: Apply the mask to the original image to get the ROI
    maskedImg = bsxfun(@times, img, cast(finalMask, 'like', img));

    % Step 6: Compute average RGB values for the ROI
    R = mean(maskedImg(:,:,1)(finalMask), 'all');
    G = mean(maskedImg(:,:,2)(finalMask), 'all');
    B = mean(maskedImg(:,:,3)(finalMask), 'all');

    % Step 7: Output average R, G, and B values to the results cell array
    results{k, 1} = barcodeID;
    results{k, 2} = R;
    results{k, 3} = G;
    results{k, 4} = B;
end

% Convert the cell array to a table
resultsTable = cell2table(results, 'VariableNames', {'barcodeID', 'R', 'G', 'B'});

% Write the table to a CSV file
writetable(resultsTable, 'matlab_leaf_rgb.csv');
