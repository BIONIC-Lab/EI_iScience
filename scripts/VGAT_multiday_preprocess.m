%code to load in VGAT5 data all recorded days
%% this script loads the calcium data for each animal and preprocesses it a bit


% the script will prompt you to select the directory for each animal
% then it will do the necessary quantification of calcium data
% run this cell to create the struct 'stimdata'


% note this script requires matlab 2019 or later
cd('C:\Users\chugh\Documents\GitHub\ExcitatoryInhibitory')


% input the animal numbers and initialize the structure
animalids = {'EIF02','EIM01','EIM09','VGAT5'};
stimdata = init_data(animalids);

%% functions

function stimdata = init_data(animalids)
        
    curdir = pwd;
    
    stimdata = struct();
    
    
    
    for animalnum = 1:length(animalids)
        
        % set animalid
        tmpid = animalids{animalnum};
        stimdata(animalnum).ID = tmpid;
        
        % save the directory
        waitfor(msgbox(sprintf('Select directory for animal number %s',tmpid)))
        stimdata(animalnum).maindir = uigetdir();
        
        %need to set the animal ID based on what was selected
        %stimdata(animalnum).ID = extractAfter(stimdata(animalnum).maindir, 'E:\ChronicEI\');
        
        % option to add animal information (electrode id, notes on how animal was recorded, surgery info etc)
        % to do this we would need to set a structure for saving
        % information about the animal note this should also include the
        % order in which the data are organized for consistency
        
        % get the roiinformation for each day

        stimdata = getroidata(stimdata,animalnum);
        
    end
    
    cd(curdir)

end

function stimdata = getroidata(stimdata,animalnum)
    
    % save the working directory
    curdir = pwd;
    
    %enter the main dir if not already there
    cd(stimdata(animalnum).maindir)
    
    % determine which days currently have information
    cd('datafiles')
    dayinfo = dir('d*');
    daydirs = {dayinfo.name};
    
    days = zeros(length(daydirs),1);
    for iday = 1:length(daydirs)
        
        days(iday) = str2num(daydirs{iday}(2:end));
        
    end
    [days,sortinfo] = sort(days);
    
    for iday = 1:length(days)
        
       stimdata(animalnum).data(iday).day = daydirs{sortinfo(iday)}; 
        
    end
    
    
    % get the trode center for each day
        % note we do not currently have trode center information so this is
        % inaccurate
        
    stimdata = gettrodecenter(stimdata,animalnum);
    
    % once trode center is input we can calculate distance
    % stimdata = calcdistance(stimdata,animalnum);
    
    
    % get the excitatorya nd inhibitory mask data
    stimdata = getEImask(stimdata,animalnum);
   
    % get the traces for each stimulation and preprocess them
    stimdata = getcalciumdata(stimdata,animalnum);
    
    % determine if cells were activated
    stimdata = getactivationmask(stimdata,animalnum);
    
    % determine if cells were inhibited/depressed
    stimdata = inhibitionmask(stimdata,animalnum,'thresh');
    
    % ss/on
    stimdata = sson(stimdata,animalnum);
    
    % activationtime
    stimdata = activationtime(stimdata,animalnum);
    

end

function stimdata = gettrodecenter(stimdata,animalnum)
    
    % we do not have the trode center for each day
        % just put it right in the middle for now
    
    for iday = 1:length(stimdata(animalnum).data)
        
        stimdata(animalnum).data(iday).trodecenter = {'not defined'};
        
    end

end

function stimdata = getEImask(stimdata,animalnum)

    curdir = pwd;
    
    for iday = 1:length(stimdata(animalnum).data)
        % enter the directory that has the information
        
        day = stimdata(animalnum).data(iday).day;        
        roiloc_dir = fullfile(stimdata(animalnum).maindir,'datafiles',day,'roifiles','locationdata');
        cd(roiloc_dir)

        % get the proper file name 
        fileinfo = dir('STD*COMEllipse.csv');
        filename = fileinfo.name;

        % read the file
        data1 = readtable(filename);

        if strcmp(stimdata(animalnum).ID, 'VGAT5')

            % karthik went the long way 'round to mark excitatory vs inhibitory
            neurons = data1{:,10};

            [imask,emask] = deal(zeros(length(neurons),1));


            eindex = find(contains(neurons,'E'));
            emask(eindex) = 1;
            imask(~emask) = 1;

        else
            % slice 1 is inhibitory image
            neurons = data1{:,8};
            [imask,emask] = deal(zeros(length(neurons),1));

            Iindex = neurons == 1;
            imask(Iindex) = 1;
            emask(~Iindex) = 1;

        end

        stimdata(animalnum).data(iday).emask = emask;
        stimdata(animalnum).data(iday).imask = imask;
    end

    cd(curdir)
        


end

function stimdata = getcalciumdata(stimdata,animalnum)

    estims = {'10u','tbs','10b','100'};
    %vstims = {'v30','v10'};
    
    for iday = 1:length(stimdata(animalnum).data)
        
       % set a field so we always know the order
       stimdata(animalnum).data(iday).Estimorder = estims;
       %stimdata(animalnum).data(iday).vstimorder = vstims;
       
       % do separate analyses for vstim and estim
       
       stimdata(animalnum).data(iday).estim.dff = loadrawcalciumtrace(stimdata,animalnum,iday,'estim');
       %stimdata(animalnum).data(iday).vstim.dff = loadrawcalciumtrace(stimdata,animalnum,iday,'vstim');
        
    end

end

function calciumdata = loadrawcalciumtrace(stimdata,animalnum,iday,stimtype)

    % this project has a 10s delay followed by 30 s estim followed by 20
    % post stim period with a frame rate of ~30hz most likely 0.034 frate
    % and we set the number of frames to be 1800
    
    % vstim was either delivered the same way, or off for 15 s on ofor 5
    
    % enter the right directory
    curdir = pwd;
    maindir = stimdata(animalnum).maindir;
    day = stimdata(animalnum).data(iday).day;
    roicalciumdata = fullfile(maindir,'datafiles',day,'roifiles', 'calciumdata');
    cd(roicalciumdata);
    
    % get the csv file names
    csvinfo = dir([stimdata(animalnum).data(iday).day '*.csv']);
    filenames = {csvinfo.name};
    
    % make not of the number of rois for initialization
    numrois = length(stimdata(animalnum).data(iday).emask);
    
    if strcmp(stimtype,'estim')
       
        %initialize data
       numstims = length(stimdata(animalnum).data(iday).Estimorder);
       numframes = 1800;
       tmpdata = zeros(numframes,numrois,numstims); % frames,number of neurons, 4 estim parameters
       
       for istim = 1:numstims
           
           % get appropriate filename
           stimname = stimdata(animalnum).data(iday).Estimorder{istim};
           fileidx = find(contains(filenames,stimname));
           filename = filenames{fileidx};
           
           % load the file and remove first column containing the frame idx
           tmpfile = csvread(filename,1,0); % changed to csvread 9/10/21 kcs
           tmpdata(:,:,istim) = tmpfile(:,2:end);
           
       end
       
       % filter the data and calculate dff
       tmpdata = filterdata(tmpdata);
       calciumdata = getdff(tmpdata,stimtype);
        
    else
        
        % v stim has two separate lengths
        tmpdata = cell(2,1);
        tmpvstim30 = zeros(1800,numrois);
        tmpvstim10 = zeros(6150,numrois);
        
        % get the 30s stim first
        fileidx = find(contains(filenames,'v30'));
        filename = filenames{fileidx};
        tmpfile = readmatrix(filename);
        tmpvstim30 = tmpfile(:,2:end);
        
        tmpvstim30 = filterdata(tmpvstim30);
        tmpdata{1} = getdff(tmpvstim30,'estim'); % keep it vstim since the delay is the same here
        
        % get the 10 trials
        fileidx = find(contains(filenames,'v10'));
        filename = filenames{fileidx};
        tmpfile = readmatrix(filename);
        tmpvstim10 = tmpfile(:,2:end);
        
        tmpvstim10 = filterdata(tmpvstim10);
        tmpdata{2} = getdff(tmpvstim10,stimtype);
        
        calciumdata = tmpdata;
        
    end
    
    cd(curdir)

end

function filtcells = filterdata(tmpdata)

    % construct the filter
    frate = 1/0.034;
    cutoff_freq = 5/frate;
    b = fir1(10,cutoff_freq,'low');
    
    % filter the data
    filtcells = single(filtfilt(b,1,tmpdata));

end

function dff = getdff(tmpdata,stimtype)

    if strcmp(stimtype,'estim')
        stimdelay = 10;
        frate = 0.034;
        baseframes = ceil(stimdelay/frate);
        
        base = mean(tmpdata(1:baseframes,:,:),1);
        base = repmat(base,size(tmpdata,1),1,1);
        
        dff = (tmpdata-base)./base;
        
        
    else
        
        stimdelay = 15;
        frate = 0.034;
        
        baseframes = ceil(stimdelay/frate);
        
        base = mean(tmpdata(1:baseframes,:,:),1);
        base = repmat(base,size(tmpdata,1),1,1);
        
        dff = (tmpdata-base)./base;
    end

end


function stimdata = getactivationmask(stimdata,animalnum)

    % do estim
        frate = 0.034;
        stimdelay = 10;
        stimduration = 30;
        
        
        
        stdmult = 3;
        
        
        for iday = 1:length(stimdata(animalnum).data)
            
           
            
            % get the data and thresh
            tmpdata = stimdata(animalnum).data(iday).estim.dff;
            
            thresh = getthresh(tmpdata,stimdelay,frate,stdmult);
            
            % is the roi above threshold for one second total minimum?
            stimframes = [ceil((stimdelay)/frate)+1:ceil((stimdelay+stimduration)/frate)];
            framesabovethresh = 30;
            stimdata(animalnum).data(iday).estim.activationmask = abovethreshcount(tmpdata,thresh,stimframes,framesabovethresh);
            
            
            % do vstim *** note this keeps the delay to 10 seconds
            % stimframes are the same for vstim30
            
%              % set up initialization for activation mask (each day) note,
%              % this will be a neurons x 2 array for the vstim
%             numneurons = length(stimdata(animalnum).data(iday).emask);
%             numstimconditions = 2; % vstim 30 and vstim10
%             tmpactivationmask = zeros(numneurons,numstimconditions);
%             
%             tmpdata = stimdata(animalnum).data(iday).vstim.dff{1};
%            
%             thresh = getthresh(tmpdata,stimdelay,frate,stdmult);
%            
%             tmpactivationmask(:,1) = abovethreshcount(tmpdata,thresh,stimframes,framesabovethresh);
%             
%             % now do the 10 trials activation
%             tmpdata = stimdata(animalnum).data(iday).vstim.dff{2};
%             thresh = getthresh(tmpdata,stimdelay,frate,stdmult);
%             tmpactivationmask(:,2) = abovethreshcount(tmpdata,thresh,stimframes,framesabovethresh);
%             
%             stimdata(animalnum).data(iday).vstim.activationmask = tmpactivationmask;
            
        end
        
    %%
    
    
    

end

function thresh = getthresh(tmpdata,stimdelay,frate,stdmult)
    %% calculate the threshold and reshape so each time point is threshed
    numframes = size(tmpdata,1);

    baseframes = ceil(stimdelay/frate);
    
    basemean = mean(tmpdata(1:baseframes,:,:),1);
    basestd = std(tmpdata(1:baseframes,:,:),[],1);
    
    thresh = basemean+stdmult*basestd;
    
    thresh = repmat(thresh,numframes,1,1);

end

function activationmask = abovethreshcount(tmpdata,thresh,activationframes,minframesabovethresh)
%% determine which neurons have activity above thresh for a full second 
    abovethresh = tmpdata>thresh;
    abovethresh30 = squeeze(sum(abovethresh(activationframes,:,:),1));
    activationmask = abovethresh30>minframesabovethresh;
    

end

function stimdata = characterize_traces(stimdata,animalnum,iday,istim)
    %%
    % initialize data
    tmpdata = stimdata(animalnum).data(iday).estim.dff(:,:,istim);
    activatedneurons = find(stimdata(animalnum).data(iday).estim.activationmask(:,istim));
    
    % initialize mask type thing and plot stuff
    responsetype = zeros(size(stimdata(animalnum).data(iday).emask));
    time = [0:1799]*0.034-10;
    
    % initialize command window prompt
    prompt = 'Please input the response profile: \n    1. Fast reduction in activity (ON) \n    2. Slow reduction in activity \n    3. Parabolic activation (increase then decrease)\n    4. Stable activation (SS) \n    5. Increasing over time\n    6. Reactivation (fast reduction and return to high activity)\n\n';
    
    % loop through all of the traces and plot them
    for ineuron = 1:length(activatedneurons)
        
       theneuron = activatedneurons(ineuron);
       neuronthresh = getthresh(tmpdata(:,theneuron),10,0.034,3);
%        ymin = mean(tmpdata(1:300,theneuron))-mean(neuronthresh); % get a mininum value for plotting
%        ymax = max(tmpdata(301:1200,theneuron));
       
       plot(time,tmpdata(:,theneuron),'LineWidth',2,'Color','r');
       hold on
       plot(time,neuronthresh','LineWidth',2,'color',[0.6 0.6 0.6],'LineStyle',':')
       axis tight
       prettyaxes('xlabel','time','ylabel','df/f','xlim', [-10 50])
       hold off
       responsetype(theneuron) = input(prompt);
        
    end
    
    stimdata(animalnum).data(iday).estim.responsetype = responsetype;

end

function stimdata = sson(stimdata, animalnum)
%% use thresholding to initially characterize which neurons are ss vs on
    stimdelay = 10;
    stimdur = 30;
    frate = 0.034;
    
    ontime = ceil([stimdelay:frate:stimdelay+2]/frate);
    sstime = ceil([stimdelay+stimdur-2:frate:stimdelay+stimdur]/frate);
    
    for iday = 1:length(stimdata(animalnum).data)
       tmpdata = stimdata(animalnum).data(iday).estim.dff;
       thresh = getthresh(tmpdata,stimdelay,frate,3);
       
       % get the mean over the last two seconds
       thresh = squeeze(mean(thresh(sstime,:,:),1));
       tmpdata = squeeze(mean(tmpdata(sstime,:,:),1));
       abovethresh_ss = tmpdata>thresh;
       
       ssmask = stimdata(animalnum).data(iday).estim.activationmask & abovethresh_ss;
       onmask = ~ssmask;
       
       stimdata(animalnum).data(iday).estim.ssmask = ssmask;
       stimdata(animalnum).data(iday).estim.onmask = onmask;
       
    end

end

function stimdata = activationtime(stimdata,animalnum)
%% how long does it take to reach half the total calcium activity
    stimdelay = 10;
    stimdur = 30;
    frate = 0.034;
    
    stimtime = ceil([stimdelay:frate:stimdelay+stimdur]/frate);
    
    for iday = 1:length(stimdata(animalnum).data)
       
        tmpactivationtime = zeros(size(stimdata(animalnum).data(iday).estim.activationmask));
        
       tmpdata = stimdata(animalnum).data(iday).estim.dff;
       cumdf = repmat(sum(tmpdata(stimtime,:,:),1)/2,length(stimtime),1,1); 
       cumsumdf = cumsum(tmpdata(stimtime,:,:),1);
       
       halfmax = cumsumdf>cumdf;
  
       for istim = 1:4
           
          for ineuron = 1:size(halfmax,2)
              
             tmpactivationtime(ineuron,istim) = find(halfmax(:,ineuron,istim),1)*frate; 
              
          end
           
       end
       
       tmpactivationtime(~stimdata(animalnum).data(iday).estim.activationmask) = nan;
       
       stimdata(animalnum).data(iday).estim.activationtime = tmpactivationtime;
       
    end

end

function thresh = negativethresh(tmpdata, stdmult)
    %% calculate the threshold and reshape so each time point is threshed
    
    % timing variables
    stimdelay = 10;
    stimdur = 30;
    frate = 0.034; % will need to change for each animal each day
    numframes = size(tmpdata,1);
    baseframes = ceil(stimdelay/frate);
    
    % mean and std
    basemean = mean(tmpdata(1:baseframes,:,:),1);
    basestd = std(tmpdata(1:baseframes,:,:),[],1);
    
    % subtract std because below baseline may mean inhibited/depressed
    thresh = basemean-stdmult*basestd;
    
    % reshape for each time point
    thresh = repmat(thresh,numframes,1,1);


end

function stimdata = inhibitionmask(stimdata,animalnum, quantmethod)
%% determine which neurons go below baseline

    % timing variables
    stimdelay = 10;
    stimdur = 30;
    frate = 0.034;
    stimtime = ceil([stimdelay:frate:stimdelay+stimdur]/frate);
    
    %loop and days
    
        
    for iday = 1:length(stimdata(animalnum).data)

        tmpdata = stimdata(animalnum).data(iday).estim.dff;
        tmpactivationmask = stimdata(animalnum).data(iday).estim.activationmask;

        if strcmp(quantmethod,'thresh')
            % do the same thresholding for activation but below baseline
            minframesbelowthresh = 30;
            thresh = negativethresh(tmpdata,3);

            belowthresh = tmpdata<thresh;
            belowthresh30 = squeeze(sum(belowthresh(stimtime,:,:),1));
            tmpinhibitionmask = belowthresh30>minframesbelowthresh;

            stimdata(animalnum).data(iday).estim.inhibitionmask = tmpinhibitionmask;

        else
            % sum calcium data below baseline -> inhibited
            threshval = 0;

            tmpsum = squeeze(sum(tmpdata(stimtime,:,:),1));
            tmpinhibitionmask = tmpsum<threshval;

            stimdata(animalnum).data(iday).estim.inhibitionmask = tmpinhibitionmask;

        end
        
    end
    
    

end
