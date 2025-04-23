function activationmask = getactivationmask(data,frate,basetime,stimtime,threshmult,activated)
% usage: activationmask = getactivationmask(data,frate,basetime,threshmult,activated)
% 
% 
% takes the data in as filter cellular data as df/f. computess the baseline mean and stantard deviation
% Returns a logical array of which cells had calcium activity that exceeded the baseline + threshmult * basestd
%    note this function should be modified for default values of the
%    various stimulation conditions so you can input data, and varargin
% inputs
%     data        - time x cells (assumes to be filtered and computed as df/f0)
%     frate       - frame period;
%     basetime    - time in seconds for the baseline period
%     stimtime    - time in seconds for the stimulaiton period
%     threshmult  - numeric how many standard deviations above the mean the activity is
%     activated   - logical to return the activated cells defaults to true,
%                   otherwise it determines whether a cell was inhibited by beling below
%                   the threshold instead of above
% 
% output
%    activationmask - logical cells by 1


    % determine how many frames to average for baseline
    baseframes = floor(basetime/frate);
    basethresh = mean(data(1:baseframes,:),1) + threshmult * std(data(1:baseframes,:),[],1);
    
    % get the stimulation frames
    stimframes = floor(basetime/frate)+1:ceil((basetime+stimtime)/frate);

    % Determine how many frames the cellular activity is above its
    % respective threshold
    if ~exist('activated','var')
        activated = true;
    end

    if activated

        thresholded = data(stimframes,:) > repmat(basethresh,length(stimframes),1);
    else

        thresholded = data(stimframes,:) < repmat(basethresh,length(stimframes),1);
    end

    % return a logical array when the data crossed threshold for at least
    % 30 frames during stimulation - note we may want to modify this
    activationmask = squeeze(sum(thresholded,1)>30);


    