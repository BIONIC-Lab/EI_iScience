function [emask,imask,exclmask,hist_temp,meangreen,meanred,opt_perc] = geteimask(basedir,animalid,rmsmult,basetime,frate,plot_func,stimpatterns,des_length)
% usage eimask = geteimask(basedir,animalid)
% 
% This function loads in the red and green roidata for all stimulation conditions
% Then computes the average value during baseline
% Fits a linear model between the red and green baselilne activity in order to correct for the green fluorescence contaminating the red channel
% corrects the red value given the fit
% 
% returns 3 logical arrays corresponding to where the average red fluorescence across all 4 trains lies 
% relative to the regression line
% Excitatory cells are any cells with red fluorescence below the regression line
% Inhibitory cells are any cells with red fluorescence above the line plus root mean square * rms mult 
% excluded cells (exclmask) are the cells that are above the regression line but below the rms*rmsmult
% 
% Inputs
% 
%         basedir         - path to the github directory for the data
%         animalid        - the animalid for going into the right folder
%         rmsmult         - integer for how many rms above the line to exclude cells
%         basetime        - time in seconds to calculate the baseline
%         frate           - the frame period
% 
% Outputs
% 
%         emask           - logical array cells x 1 where any cell with red fluorescence below regression line
%         imask           - logical array cells x 1 where any cell with red fluorescece above the regression line + rms mult*rms
%         exclmask        - any remaining cell that didnt fit the above criteria
% 
% 


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
        
        %make desired length
        tmpred =  [tmpred; NaN(des_length-size(tmpred,1),size(tmpred,2))];
        tmpgreen =  [tmpgreen; NaN(des_length-size(tmpgreen,1),size(tmpgreen,2))];

        redvals = cat(3,redvals,tmpred(:,2:end)); %this won't work if the ROI maps are not consistent
        greenvals = cat(3,greenvals,tmpgreen(:,2:end));

    end

    % average across stim patterns and baseline
    onset = floor(basetime/frate);
    baseframes = 1:onset;
    stimframes = onset+1:onset+floor(3/frate); %use the first 3 s
    
    
    meanred = squeeze(mean(mean(redvals(baseframes,:,:),3),1));
    std_red = squeeze(std(mean(redvals(baseframes,:,:),3),1));
    meangreen = squeeze(mean(mean(greenvals(baseframes,:,:),3),1));
    
    greenvals_norm = (greenvals(:,:,1)-meangreen)./meangreen;
    
%     %just focus on 100 Hz
%     meanred = squeeze(mean(redvals(baseframes,:,3),1));
%     meangreen = squeeze(mean(greenvals(baseframes,:,3),1));
%     
%     %can use diff to find a ratio
%     activered = squeeze(mean(redvals(stimframes,:,3),1));
%     activegreen = squeeze(mean(greenvals(stimframes,:,3),1));
%     
%     %difference between first 3 s of evoked activity and baseline
%     diffred = activered-meanred;
%     diffgreen = activegreen-meangreen;
%     
%     %percent of contamination i.e. what percent of increase in green
%     %channel makes it to red channel
%     perc_cont = diffred./diffgreen; 
%     perc_cont(perc_cont>2) = 2; %assume anything beyond this is an outlier based on small evoked activity
%     perc_cont(perc_cont<0) = 0;  %green channel should never negatively influence red channel
    %out_idx = isoutlier(perc_cont);
    %perc_cont(out_idx) = []; %remove big outliers for easier histogram visualization
    
    %can now try to calculate the intensity of green that makes it into red
    %unfortunately this doesn't seem to work, green and red channel have
    %stronger influence on each other than expected from this
%     meanred_fix = meanred-(meangreen.*perc_cont);
%     order_idx = sort(meanred_fix);
%     inh_thres = order_idx(ceil(length(meanred_fix)*0.8));
%     exc_thres = order_idx(ceil(length(meanred_fix)*0.6));
%     imask = meanred_fix > inh_thres;
%     emask = meanred_fix < exc_thres;

    % create a linear model between the red and the green values
    mdl = fitlm(meangreen,meanred);
    intercept = table2array(mdl.Coefficients(1,1));
    slope = table2array(mdl.Coefficients(2,1));

%     intercept = 0;
%     slope = table2array(mdl.Coefficients(1,1));
    %% method for drop all points below regression line, recalculate "Method 1"
%     meanred_drop = meanred';
%     meanred_drop(meanred'<mdl.Fitted) = NaN;
%     meangreen_drop = meangreen';
%     meangreen_drop(meanred'<mdl.Fitted) = NaN;
%     
%     mdl_drop = fitlm(meangreen_drop,meanred_drop);
%     imask = meanred' > mdl_drop.Fitted;
%     emask = meanred' < mdl.Fitted;
    
%% save for histogram generation across animals 
    hist_temp = meanred./meangreen;
%     intercept = table2array(mdl.Coefficients(1,1))/mean(meangreen);
%     
%     %use the intercept from the regression
%     
%     %just use red to green ratio and divide into percents within each
%     %animal
%     hist_order = sort(hist_temp,'descend');
%     idx = ceil(length(hist_order)*0.2);
%     cut_off = hist_order(idx);
%     idx2 = ceil(length(hist_order)*0.6);
%     cut_off2 = hist_order(idx2);
    
%% method for dividing by percent expected and adjusting slope
%     diff_area = 0;
%     for perc = 0.05:0.01:0.5 %finding 'optimal' percent
%         slope = 1.01;
%         imask = meanred' > mdl.Fitted;
%         while sum(imask) >= perc*length(meanred)
%             imask = meanred' > meangreen'*slope+intercept;%mdl.Fitted*slope;
%             slope = slope+0.1;
%         end
%         emask = meanred' < mdl.Fitted;
%         diff_temp = nanmean(nanmean(greenvals_norm(:,imask),2) - nanmean(greenvals_norm(:,emask),2)); %only using 100 Hz for now, can switch between 
%         if diff_temp > diff_area
%             opt_perc = perc;
%             diff_area = diff_temp;
%             imask_keep = imask;
%             emask_keep = emask;
%             slope_keep = slope;
%         end
%     end
% 
%     imask = imask_keep;
%     emask = emask_keep;
%% method for dividing by percent expected and adjusting intercept
%     diff_area = 0;
%     for perc = 0.05:0.01:0.5 %finding 'optimal' percent
%         intercept = 0;
%         imask = meanred' > mdl.Fitted;
%         while sum(imask) >= perc*length(meanred)
%             imask = meanred' > meangreen'*slope+intercept;%mdl.Fitted*slope;
%             intercept = intercept+1;
%         end
%         emask = meanred' < mdl.Fitted;
%         diff_temp = nanmean(nanmean(greenvals(:,imask,3),2) - nanmean(greenvals(:,emask,3),2)); %only using 100 Hz for now, can switch between      
%     end
    
%% method replotting data on origin and finding standard deviation "Method 3"
%     mdl_origin = meanred'-mdl.Fitted; %corrected red intensity for model
%     std_origin = nanstd(mdl_origin);
%     mean_origin = nanmean(mdl_origin);
% %     
%     %if we want to cut-off at standard deviation (will make threshold
%     %around 16%)
%     imask = mdl_origin'>mean_origin+std_origin;
%     emask = mdl_origin'<mean_origin;
%     
%% can also choose data points in the model so that we get 20% inhibitory "Method 2"
%     mdl_origin = meanred'-mdl.Fitted; %corrected red intensity for model
%     mean_origin = nanmean(mdl_origin);
%     n=numel(meangreen);
%     xy3=[meangreen' meanred' mdl_origin];
%     y4=sortrows(xy3,3);
%     thresh2=ceil(round(n*0.90)); %upper 20%
% %     ii=0;
% %     for i=thresh2:n
% %         ii=ii+1;
% %         x5(ii)=y4(i,1);
% %         y5(ii)=y4(i,2);
% %     end
%     cutoff = y4(thresh2,3);
%     imask = mdl_origin'>cutoff;
%     emask = mdl_origin'<mean_origin;

%% alt method using fraction of red to green
%     imask = hist_temp >= cut_off+intercept;
%     emask = hist_temp <= cut_off2+intercept; %this was approximately the value I found

%% original method
%     %add in additional criteria
%     meanred_sort = sort(meanred);
%     threshold = meanred_sort(ceil(length(meanred_sort)*0.8));
%     
%     meangreen_sort = sort(meangreen);
%     threshold2 = meangreen_sort(ceil(length(meangreen_sort)*0.4));

    imask = meanred' > mdl.Fitted + rmsmult*mdl.RMSE; % & meanred'>threshold & meangreen'>threshold2;
    emask = meanred' < mdl.Fitted; % & meanred'<threshold & meangreen'>threshold2;

    %find upper percent based on red channel
%     meanred_sort = sort(meanred);
%     length_red = length(meanred);
%     top_perc = ceil(length_red*0.9);
%     threshold = meanred_sort(top_perc);
%     bot_perc =  ceil(length_red*0.4);
%     threshold2 = meanred_sort(bot_perc);
% 
% 
%     imask = meanred' >= threshold;
%     emask = meanred' <= threshold2;

    exclmask = ~imask & ~emask;
    
    if plot_func
        %plot histogram before correction
        figure
        h1 = histogram(meanred,50);
        xlabel('Baseline F Red Channel')
        ylabel('# of ROIs')
        title('Red Channel ROIs before correction')
        
        %plot cumulative sum and divisions
        yyaxis right
        cs = cumsum(h1.Values);
        plot(h1.BinEdges(2:end),cs, 'LineWidth', 2)
        ylabel('Cumulative sum of ROIs')

        %plot model
        figure
        
        %mark exc, inh, and mix to color     
        scatter(meangreen(emask), meanred(emask), 'filled', 'b')
        hold on
        scatter(meangreen(imask), meanred(imask), 'filled', 'r')
        scatter(meangreen(exclmask), meanred(exclmask), 'filled', 'm')

%         scatter(meangreen(emask), mdl_origin(emask), 'filled', 'b')
%         hold on
%         scatter(meangreen(imask), mdl_origin(imask), 'filled', 'r')
%         scatter(meangreen(exclmask), mdl_origin(exclmask), 'filled', 'm')
        
        plot(meangreen, mdl.Fitted, 'k')
        plot(meangreen, mdl.Fitted+mdl.RMSE, 'c')
%         plot(meangreen, meangreen'*slope_keep+intercept, 'c')
%         plot(meangreen, mdl_drop.Fitted, 'c')
%         yline(mean_origin,'k--')
%         yline(mean_origin+std_origin,'c--')
%         yline(cutoff,'c--')
        
        legend('off')
        xlabel('Baseline F Green Channel')
        ylabel('Baseline F Red Channel')
        title('Influence of GCaMP on Red Channel')
    end


