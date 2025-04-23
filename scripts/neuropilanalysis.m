function stimdata = neuropilanalysis(extensionwidth,binsize,numframes)
%%% usage stimdata = neuropilanalysis()
%
%currently this function is in testing
gcp;

% test on one animal
animals = {'EIM08','EIM09','VGAT2','VGAT5'};
train_types = {'10 Hz', '10 Hz Burst', '100 Hz','TBS'};
trainfoldersuffix = [2940,2938,2941,2939;...
    3048,3047,3046,nan;...
    2173,2172,2174,2175;...
    2445,2442,2444,2443];


extensionMults = [1 2 3]; % multiples of the diameter of the neuron
binsize = 20; %um


% prompt user to select the local data location for loading raw images
waitfor(msgbox('Please select the base directory for day 0 processing of animal data'))
localdataloc = uigetdir();

% prompt user to select the local location for accesing the github
% repository
waitfor(msgbox("Please select the local location for the githubrepository for the EI project"))
basedir = uigetdir();
basetime = 10;
numframes = 1800; % numframes as an input parameter because vstim would be a different size

% select a savedirectory for saving the data
waitfor(msgbox("Please select the local location to save the neuropildata"))
savedir= uigetdir();

neuropildata = struct();

% load in the location information
locationdir = fullfile(basedir,'Locations');
for animalnum = 1:length(animals)
    
    neuropildata(animalnum).id = animals{animalnum};
    neuropildata(animalnum).dataorder = train_types;
    neuropildata(animalnum).stimdatadirnum = trainfoldersuffix(animalnum,:);
    neuropildata(animalnum).npilbinsize = binsize;
    neuropildata(animalnum).annuluswidth = extensionmults;


    % get the filename for the animal
    locationinfo = dir(fullfile(locationdir,sprintf('%s*Locations.csv',animals{animalnum})));
    locationinfo = readmatrix(locationinfo.name);
    locationinfo = locationinfo(2:end,2:end);
    numrois = size(locationinfo,1);

    % load in the elecrode location for the neuropil bin calculation
    trodeinfo = dir(fullfile(locationdir,sprintf('%s*Electrode*.csv',animals{animalnum})));
    trodeinfo = readmatrix(trodeinfo.name);
    trodeinfo = trodeinfo(2,2:3);
    neuropildata(animalnum).trodeloc = trodeinfo;



    % initialize some data
    [animalroidata] = deal(zeros(numframes,numrois,length(train_types)));
    animalnpilannulus = zeros(numframes,numrois,length(extensionMults),length(train_types));
    
    animalnpilbins = []; % empty becasue we dont know the number of bins and the allocation likely wont be too much
    [animalannuluscorr,animalbincorr] = deal(zeros(numrois,length(train_types)));


    for curr_train = 1:length(train_types)
        
        if ~isnan(trainfoldersuffix(animalnum,curr_train))
        % load in the green channel images 
        imageinfo = dir(fullfile(localdataloc,animals{animalnum},'datafiles','d000','2p', sprintf('Tseries*%0.f',trainfoldersuffix(animalnum,curr_train))));
        [tmpseries,~,imgdata] = doimport2p(fullfile(imageinfo.folder,imageinfo.name),3);
        [h,w,t] = size(tmpseries);


        % get the indices of the neuropil
        %extensionwidth_pxl = floor(extensionwidth/imgdata.xyres);
        [roiindices,annulusindices,annulusmask] = getneuropilannulus(std(tmpseries,[],3),locationinfo,extensionwidth_pxl,false);
        

        % extend linear roi indices to include all frames - faster but
        % identical to creating a single frame nask and multiplying it by
        % all frames
        roiindices_ext = cell(size(roiindices));
        annulusindices_ext = cell(size(annulusindices));
        for iroi = 1:numrois
            roiindices_ext{iroi} = bsxfun(@plus,roiindices{iroi},(0:t-1)*w*h);
            for idilation = 1:length(extensionMults)
                annulusindices_ext{iroi,idilation} = bsxfun(@plus,annulusindices{iroi,idilation},(0:t-1)*w*h);
            end
        end

        [rawcelldata] = deal(zeros(t,numrois));
        rawneuropildata = zeros(t,numrois,length(extensionMults));
        for iroi = 1:numrois

            roivalues = tmpseries(roiindices_ext{iroi});
            
            rawcelldata(:,iroi) = mean(roivalues,1);

            for idilation = 1:length(extensionMults)
                neuropilvalues = tmpseries(annulusindices_ext{iroi,idilation});
                rawneuropildata(:,iroi,idilation) = mean(neuropilvalues,1);
            end

        end

        % filter the cellular and neuropil data - smooth it as chris did
        % but maybe construct the filter later
        filtcells = getdff(smoothdata(rawcelldata,'Gaussian',[7,7]),basetime,imgdata.frate);
        filtneuropil = zeros(numrois,length(extensionMults));
        for idilation = 1:length(extensionMults)

            filtneuropil(:,idilation) = getdff(smoothdata(rawneuropildata(:,:,idilation),'Gaussian',[7,7]),basetime,imgdata.frate);

        end
        

        % calculate the correlation between each aroi and its corresponding
        annuluscorr = zeros(numrois,length(extensionMults));
        for idilation = 1:length(extensionMults)

            annuluscorr(:,idilation) = diag(corr(filtcells,filtneuropil));

        end

        %=================================================================
        
        % use roi location/distance data loaded above to calculate the
        % correlation with the respective radial bin
        
        % calculate the distance from the electrode in mucrons
        roidistances = pdist2(locationinfo(:,1:2),trodeinfo,'euclidean').*imgdata.xyres;
        neuropildata(animalnum).roidist = roidistances;
        
        
        % create a distance map of distances from the electrode
        troderoi = images.roi.Circle([],'Center',trodeinfo,'Radius',30/2/imgdata.xyres);
        trodemask = createMask(troderoi,h,w);
        distmask = bwdist(trodemask.*imgdata.xyres);
        
        maxdist = max(distmask(:))+binsize;

        % get the bin indices and calculate the neuropil activity in the
        % bins
        bins = [0:binsize:maxdist];
        neuropildata(animalnum).bins = bins;
        neuropilbins = zeros(t,length(bins) - 1);
        for ibin = 1:length(bins) - 1

            tmpimg = false(size(distmask));
            tmpimg(distmask>bins(ibin) & distmask<=bins(ibin+1)) = true;
            
            tmpbinindices = bsxfun(@plus,find(tmpimg),(0:t-1)*w*h);
            
            neuropilbins(:,ibin) = mean(tmpseries(tmpbinindices));
        
        end

        % filter the neuropil data
        neuropilbins = getdff(smoothdata(neuropilbins,'gaussian'),basetime,imgdata.frate);

        % determine which bins each roi is in
        [~,~,roibinloc] = histcounts(roidistances,bins);

        % calculate the correlation bewteen all rois in each bin
        bincorr = zeros(size(annuluscorr));
        for ibin = 1:max(roibinloc)

            tmproiindex = roibinloc == ibin;
            if any(tmproiindex)
                bincorr(tmproiindex) = corr(filtcells(:,tmproiindex),neuropilbins(:,ibin));
            end

        end

        % populate the full arrays with the relevant data
        animalroidata(:,:,curr_train) = filtcells;
        animalnpilannulus(:,:,:,curr_train) = filtneuropil;
        animalnpilbins = cat(3,animalnpilbins,neuropilbins); 

        animalannuluscorr(:,curr_train) = annuluscorr;
        animalbincorr(:,curr_train) = bincorr;
    
        neuropildata(animalnum).frate = imgdata.frate;
        neuropildata(animalnum).xyres = imgdata.xyres;
        end
    end
    
    
    

    neuropildata(animalnum).cellcalcium = animalroidata;
    neuropildata(animalnum).npilannulus = animalnpilannulus;
    neuropildata(animalnum).npilbins = animalnpilbins;

    neuropildata(animalnum).binedges = bins;

    neuropildata(animalnum).annuluscorr = animalannuluscorr;
    neuropildata(animalnum).bincorr = animalbincorr;


end

cd(savedir)
save('neuropildata.mat','neuropildata');
