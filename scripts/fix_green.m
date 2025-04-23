%test script to correct green channel instead of red channel
function [mdl_fit,intercept] = fix_green(basedir,animalid,basetime,frate,stimpatterns)


    % navigate to the right animal directory and load all of the roi files
    % for green and red
    
    opt_perc = 0; %default

    redvals = [];
    greenvals = [];

    for istim = 1:length(stimpatterns) %use all trains for consistency

        redinfo = dir(fullfile(basedir,'Averages',sprintf('%s',stimpatterns{istim}),'Red',sprintf('%s*.csv',animalid)));
        greeninfo = dir(fullfile(basedir,'Averages',sprintf('%s',stimpatterns{istim}),'Green',sprintf('%s*.csv',animalid)));

        tmpred = readmatrix(fullfile(redinfo.folder,redinfo.name));
        tmpgreen = readmatrix(fullfile(greeninfo.folder,greeninfo.name));

        redvals = cat(3,redvals,tmpred(:,2:end)); %this won't work if the ROI maps are not consistent
        greenvals = cat(3,greenvals,tmpgreen(:,2:end));

    end

    % average across stim patterns and baseline
    onset = floor(basetime/frate);
    baseframes = 1:onset;
    stimframes = onset+1:onset+floor(3/frate); %use the first 3 s   
    
    meanred = squeeze(mean(mean(redvals(baseframes,:,:),3),1));
    meangreen = squeeze(mean(mean(greenvals(baseframes,:,:),3),1));    

    % create a linear model between the red and the green values
    mdl = fitlm(meanred,meangreen);
    
    mdl_fit = mdl.Fitted;
    intercept = table2array(mdl.Coefficients(1,1));
