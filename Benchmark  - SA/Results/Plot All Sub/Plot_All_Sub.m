close
clear
clc
%%
% Create a new figure
fig = figure;

% Define the size of the figure
fig.Position = [250 0 1200 1200]; % [left bottom width height]

% Define the list of filenames
fileNames = {'Results_SA_Area_try_1.png', 'Results_SA_Area_try_2_BAD.png', 'Results_SA_Area_try_3_BAD.png', 'Results_SA_Area_try_4.png', 'Results_SA_Area_try_5.png'};

for i = 1:numel(fileNames)
    subplot(3,2,i); % Create a new subplot within the figure
    img = imread(fileNames{i}); % Read the .png file
    imshow(img, 'InitialMagnification', 'fit'); % Display the image in the subplot
end
