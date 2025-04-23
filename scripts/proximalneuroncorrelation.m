function neuroncorr = proximalneuroncorrelation(neuronums,activated                 )
% usage neuroncorr = proximalneuroncorrelation(neuronnums)
% 
% this function prompts the user to select the local github directory then
% loads in the cellular location info, and the calcium data for all of the trials of all animals
% uses another custom function to determine which cells are activated, and which are excitatory or inhibitory
% then for each cell calculates its correlation with the neuronnums closest cells
% 
% then it returns a struct describing the correlation between:
% E cells and other E cells
% E cells and other I cells 
% I cells and other I cells
% 
% inputs
% 
%     neuronnums    - integer for the number of closest neurons to correlate with
% 
% outputs
% 
%     neuroncorr    - struct with fields EE EI II


    animals = {'EIM08','EIM09','VGAT2','VGAT5'};
    train_types = {'10 Hz', '10 Hz Burst', '100 Hz', 'TBS'};

    waitfor(msgbox("Please select the local location for the githubrepository for the EI project"))
    basedir = uigetdir();
    basetime = 10;
    numframes = 1800; % numframes as an input parameter because vstim would be a different size
    rmsmult = 1;
    frate = 0.034;
    stdmult = 3;


    for animalnum = 1:length(animals)
        
        % load in the cellular location information
         % get the filename for the animal
        locationinfo = dir(fullfile(locationdir,sprintf('%s*Locations.csv',animals{animalnum})));
        locationinfo = readmatrix(locationinfo.name);
        locationinfo = locationinfo(2:end,2:end);
        numrois = size(locationinfo,1);

        % get the e and i mask
        [emask,imask,other] = geteimask(basedir,animals{animalnum},rmsmult,basetime,frate);

        for curr_train = 1:length(traintypes)

            % load in the calcium data 
            calciuminfo = dir(fullfile(basedir,'Averages',train_types{curr_train},sprintf('%s*.csv',animals{animalnum})));
            tmpcalcium = readmatrix(calciuminfo.name);
            tmpcalcium = tmpcalcium(2:end,2:end);

            % filter and calculate dff
            dff = getdff(smooth(tmpcalcium,'Gaussian',[7,7]),basetime,frate);
            clear tmpcalcium

            % get the activationmask
            activationmask = getactivationmask(data,frate,basetime,threshmult,activated);

            % calculate correlation between e neurons and e neurons
            eindex = find(emask&activationmask);
            iindex = imask&activationmask;

            

            % calculate correlation between e and i neurons

            % calculate the correlation between i and i neurons

        end

    end