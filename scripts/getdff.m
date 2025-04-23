function dff = getdff(data,basetime,frate)
% usage dff = getdff(data,basetime,frate)
% 
% this function takes in data (time x cells/bins/rois) calculates the mean baseline activity
% then converts the data into dff (data - baseline)./baseline
% 
% inputs
% 
%         data        - time x roi array
%         basetime    - seconds for baseline calculation
%         frate       - frame period for collecing the images
% 
%  outputs
% 
%         dff         - time x roi array that is now in df/f instead of raw fluorescence
%         


    % calculate the mean baseline activity
    baseframes = 1:floor(basetime/frate);
    meanbase = mean(data(baseframes,:),1);

    % calculate the dff
    dff = (data - meanbase)./meanbase;