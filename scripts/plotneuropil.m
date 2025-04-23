function plotneuropil(neuropildata)

    %% initialize some things for processing
    numanimals = length(neuropildata);
    onset_start = 294;
    frate = 0.034;
    rmsmult = 1;
    stimpatterns = neuropildata(1).dataorder;
    numframes = 1800;
    basetime = 10;
    plot_func = false;

    des_length = numframes;
    neuropildata = getfijicalcium(neuropildata);
    % get the github main directory 
    waitfor(msgbox('Select the directory for the excitatoty inhibitory repository'))
    basedir = uigetdir();

    % get which cells are considered excitatory inhibitory or "mixed"
%%
    for animalnum = 1:numanimals


        tmpid = neuropildata(animalnum).id;
        % get the excitatory and inhibitory mask for each animal
        [exc_idx,inh_idx,oth_idx,hist_temp,meangreen,meanred,~] = geteimask(basedir,tmpid,rmsmult,basetime,frate,plot_func,stimpatterns,des_length);
        
        % get activationmask for refining the activation
        tmpdata = neuropildata(animalnum).fijicalcium;
        threshmult = 5;
        stimtime = 30;

        tmpactivationmask = geteiactivationmask(tmpdata,basetime,stimtime,frate,threshmult);
        
        neuropildata(animalnum).emask = exc_idx;
        neuropildata(animalnum).imask = inh_idx;
        neuropildata(animalnum).mmask = oth_idx;
        neuropildata(animalnum).activationmask = tmpactivationmask;

        % double check that correlation patterns remains the same when only
        % considering correlation during stim
        [tmpannuluscorr,tmpannulusbase] = deal(zeros(size(tmpactivationmask,1),length(neuropildata(animalnum).annuluswidth),4));
        
        for istim = 1:4

            for idilation = 1:4
                tmpneuron = neuropildata(animalnum).cellcalcium(onset_start:end,:,istim);
                tmpnpil = neuropildata(animalnum).npilannulus(onset_start:end,:,idilation,istim);
                tmpannuluscorr(:,idilation,istim) = diag(corr(tmpneuron,tmpnpil));

                tmpneuronbase = neuropildata(animalnum).cellcalcium(1:onset_start-1,:,istim);
                tmpnpilbase = neuropildata(animalnum).npilannulus(1:onset_start-1,:,idilation,istim);

                tmpannulusbase(:,idilation,istim) = diag(corr(tmpneuronbase,tmpnpilbase));
            end

        end
        neuropildata(animalnum).annuluscorrorig = neuropildata(animalnum).annuluscorr;
        neuropildata(animalnum).annuluscorr = tmpannuluscorr;
        neuropildata(animalnum).annuluscorrbase = tmpannulusbase;

    end

    % plot how the diameter of the annulus affects correlation
    % diameter increases then stays the same


    %% plot mean of cell types
    meancorrelation(neuropildata)
    %pbaspect([0.5 4.9])
    % stim p= 0.003

    %% does it change over distance
    for istim = 1:4

        plotcoor_v_dist(neuropildata,istim,1);

    end


    % all distance p<0.0001

    % do response class across cell types for each stimulation pattern

    %% how does it compare regardless of stim condition

    comparisons = plotcoor_v_dist_allstim(neuropildata,1);

    %cell type p=0.02
    %distance = 9e-34
    % no sig with all comparisions 
    % ei diff at 200+ =0.0044 with only one comparison


    %% what about baseline correlations

    comparisons = plotbasecoor_v_dist_allstim(neuropildata,1);
%%
    for istim = 1:4

        neuropil_v_responseclass(neuropildata,istim)

    end



return

function representativetrace(neuropildata,animalnum)
%%
    istim = 3;
    closethresh = 200;
    farthresh = 200;

    closeneurons = neuropildata(animalnum).roidist<=closethresh;
    farneurons = neuropildata(animalnum).roidist>=farthresh;

    emask_close = neuropildata(animalnum).activationmask(:,istim) & neuropildata(animalnum).emask & closeneurons;
    emask_far = neuropildata(animalnum).activationmask(:,istim) & neuropildata(animalnum).emask & farneurons;

    imask_close = neuropildata(animalnum).activationmask(:,istim) & neuropildata(animalnum).imask & closeneurons;
    imask_far = neuropildata(animalnum).activationmask(:,istim) & neuropildata(animalnum).imask & farneurons;
    
    tmpneurons = find(imask_far);
    for ineuron = 1:length(tmpneurons)

        plot(neuropildata(animalnum).cellcalcium(:,tmpneurons(ineuron),istim),'color','k');
        hold on
        plot(neuropildata(animalnum).npilannulus(:,tmpneurons(ineuron),1,istim),'color','r','LineStyle','--');
        title(sprintf('%0.f %.2f',tmpneurons(ineuron),neuropildata(animalnum).annuluscorr(tmpneurons(ineuron),1,istim)))
        pause(2)
        clf

    end
    %%
    eneuron_close = [92,93,94];
    eneuron_far = 107;
    ineuron_close =56;
    ineuron_far = [29,10,44];

    colors = [90,156,254;231 118 129; 157,126,213];
    colors = colors/255;

    time = [0:1799]*0.034 - 10;
    figure,
    plot(time,neuropildata(animalnum).cellcalcium(:,eneuron_close(1),istim),'LineWidth',3,'Color',colors(1,:))
    hold on
    plot(time,neuropildata(animalnum).npilannulus(:,eneuron_close(1),istim),'LineWidth',3,'Color',colors(1,:),'LineStyle','--')
    %%prettyaxes('xlabel','Time,s','ylabel','Df/f','ylim',[-0.25 3])
    legend('Neuron ROI','Neuropil')
    title(sprintf('Exc: %.3f',neuropildata(animalnum).annuluscorr(eneuron_close(1),1,istim)));
    %pbaspect([3 1 1])

%corr(neuropildata(animalnum).cellcalcium(:,eneuron_close(1),istim),neuropildata(animalnum).npilannulus(:,eneuron_close(1),1,istim))

    figure,
    plot(time,neuropildata(animalnum).cellcalcium(:,eneuron_far,istim),'LineWidth',3,'Color',colors(1,:) - 90/255)
    hold on
    plot(time,neuropildata(animalnum).npilannulus(:,eneuron_far,istim),'LineWidth',3,'Color',colors(1,:) - 90/255,'LineStyle','--')
    %%prettyaxes('xlabel','Time,s','ylabel','Df/f','ylim',[-0.25 3])
    legend('Neuron ROI','Neuropil')
    title(sprintf('Exc: %.3f',neuropildata(animalnum).annuluscorr(eneuron_far,1,istim)));
    %pbaspect([3 1 1])

    figure,
    plot(time,neuropildata(animalnum).cellcalcium(:,ineuron_close,istim),'LineWidth',3,'Color',colors(2,:))
    hold on
    plot(time,neuropildata(animalnum).npilannulus(:,ineuron_close(1),istim),'LineWidth',3,'Color',colors(2,:),'LineStyle','--')
    %%prettyaxes('xlabel','Time,s','ylabel','Df/f','ylim',[-0.25 3])
    legend('Neuron ROI','Neuropil')
    title(sprintf('Inh: %.3f',neuropildata(animalnum).annuluscorr(ineuron_close(1),1,istim)));
    %pbaspect([3 1 1])

    figure,
    plot(time,neuropildata(animalnum).cellcalcium(:,ineuron_far(1),istim),'LineWidth',3,'Color',colors(2,:) - 118/255)
    hold on
    plot(time,neuropildata(animalnum).npilannulus(:,ineuron_far(1),istim),'LineWidth',3,'Color',colors(2,:) - 118/255,'LineStyle','--')
    %%prettyaxes('xlabel','Time,s','ylabel','Df/f','ylim',[-0.25 3])
    legend('Neuron ROI','Neuropil')
    title(sprintf('Inh: %.3f',neuropildata(animalnum).annuluscorr(ineuron_far(1),1,istim)));
    %pbaspect([3 1 1])
    


return

function neuropilsubtraction(neuropildata,alpharange)

    % plot mean response of neurons at each stimulation frequency for
    % neuropul subtraction for each cell type
%%
    raw = cell(3,4);
    subtracted = cell(3,4);

    for animalnum = 1:4
        
        
        for istim = 1:4
            tmpactivationmask = neuropildata(animalnum).activationmask(:,istim);

            tmpemask = neuropildata(animalnum).emask;
            tmpimask = neuropildata(animalnum).imask;
            tmpmmask = neuropildata(animalnum).mmask;

            tmpannulus = neuropildata(animalnum).npilannulus(:,:,1,istim);
            tmpdata = neuropildata(animalnum).fijicalcium(:,:,istim);
            tmpdata(:,~tmpactivationmask) = nan;
            
            raw{1,istim} = cat(2,raw{1,istim},tmpdata(:,tmpemask));
            raw{2,istim} = cat(2,raw{2,istim},tmpdata(:,tmpimask));
            raw{3,istim} = cat(2,raw{3,istim},tmpdata(:,tmpmmask));


            tmpsubtracted_e = zeros(1800,sum(tmpemask),length(alpharange));
            for ialpha = 1:length(alpharange)

                tmpsubtracted_e(:,:,ialpha) = tmpdata(:,tmpemask) - alpharange(ialpha) .* tmpannulus(:,tmpemask); 

            end
            subtracted{1,istim} = cat(2,subtracted{1,istim},tmpsubtracted_e);

            tmpsubtracted_i = zeros(1800,sum(tmpimask),length(alpharange));
            for ialpha = 1:length(alpharange)

                tmpsubtracted_i(:,:,ialpha) = tmpdata(:,tmpimask) - alpharange(ialpha)  .* tmpannulus(:,tmpimask); 

            end
            subtracted{2,istim} = cat(2,subtracted{2,istim},tmpsubtracted_i);

            tmpsubtracted_m = zeros(1800,sum(tmpmmask),length(alpharange));
            for ialpha = 1:length(alpharange)

                tmpsubtracted_m(:,:,ialpha) = tmpdata(:,tmpmmask) - alpharange(ialpha)  .* tmpannulus(:,tmpmmask); 

            end
            subtracted{3,istim} = cat(2,subtracted{3,istim},tmpsubtracted_m);

        end

    end


    % get the labels for the activity
    tmpalldata = cell(4,3,5);

    for ialpha = 1:5
        for istim = 1:4
    
    
            for icell = 1:3
    
                tmpalldata{istim,icell,ialpha} = subtracted{icell,istim}(:,:,ialpha)';
    
            end
    
        end
    end

    tmperrorall = cellfun(@(x) std(x(:,1:294),[],2,'omitnan')',tmpalldata,'UniformOutput',false);


    for ialpha = 1:5
        ialpha = 2;
        labels = labels_during_stim(tmpalldata(:,:,ialpha),tmperrorall(:,:,ialpha));
        label_types = {'RA','SA','NA','Fac','Decrease'};
        plot_by_label(tmpalldata(:,:,ialpha),labels,label_types)

    end

    meandata = cellfun(@(x) squeeze(mean(x,2,'omitnan')),tmpalldata,'UniformOutput',false);
    % stims = neuropildata(1).dataorder;
    % neurons = {'Exc','Inh','Uncerain'};
    % for ineuron = 1:3
    % 
    %     figure,
    % 
    %     for istim = 1:4
    % 
    %         nexttile
    %         tmpraw = mean(raw{ineuron,istim},2);
    %         plot(tmpraw,'LineWidth',2,'color','k')
    %         hold on
    %         for ialpha = 1:length(alpharange)
    % 
    % 
    %             tmpsubtracted_plot = mean(subtracted{ineuron,istim}(:,:,ialpha),2);
    %             plot(tmpsubtracted_plot)
    % 
    % 
    %         end
    % 
    %         legend({'0','0.1','0.2','0.8','1'})
    %         title(sprintf('%s: %s',neurons{ineuron},stims{istim}))
    % 
    %     end
    % 
    % end


return


function neuropildata = getfijicalcium(neuropildata)
    %maindir = '/Users/kes258/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Github/ExcitatoryInhibitory/Averages';
    maindir = 'C:\Users\chugh\OneDrive\Documents\GitHub\ExcitatoryInhibitory\Averages';
    stimpatterns = neuropildata(1).dataorder;
    fulldata = cell(4,1);
    onset_start = 294; %approximate frame that ICMS started
    for istim = 1:4
        cd(fullfile(maindir,stimpatterns{istim},'Green'));
        for animalnum = 1:4
    
            tmpid = neuropildata(animalnum).id;
            fileinfo = dir(sprintf('%s*.csv',tmpid));
            tmpdata = readmatrix(fileinfo.name,'NumHeaderLines',1);

            fulldata{animalnum} = cat(3,fulldata{animalnum},tmpdata(:,2:end));
    
        end

    end

    for animalnum = 1:4
        
        tmpdata = smoothdata(fulldata{animalnum}, 'gaussian',[7,7]);
        base = mean(tmpdata(1:onset_start,:,:),1);

        neuropildata(animalnum).fijicalcium = (tmpdata - base)./base;

    end

return

function activationmask = geteiactivationmask(tmpdata,basetime,stimtime,frate,threshmult)

    baseframes = floor(basetime/frate);
    stimframes = ceil((basetime+stimtime)/frate);
    datasize = size(tmpdata);
    base = mean(tmpdata(1:baseframes,:),1);
    basestd = std(tmpdata(1:baseframes,:),[],1);

    activationmask = reshape((max(abs(tmpdata(baseframes:stimframes,:)),[],1))> (abs(base) + threshmult*basestd),datasize(2:end));


return

function neuropil_v_responseclass_v_stim(neuropildata)
    %%
    fullcorr = cell(4,3,4); % response class e/i/m istim
    [my_annulus,responselabels,my_dist] = geteilabels(neuropildata);
%%
    fulldist = cell(4,3,4);
    for iclass = 1:4


        for ineuron = 1:3


            for istim = 1:4
                
                tmpdata = my_annulus{istim,ineuron};
                tmpdist = my_dist{istim,ineuron};
                fullcorr{iclass,ineuron,istim} = tmpdata(responselabels{istim,ineuron} == iclass);
                fulldist{iclass,ineuron,istim} = tmpdist(responselabels{istim,ineuron} == iclass);
            end

        end

    end


    stimcolors_dark = [146,35,90; 172, 86, 14;27 114 47 ;217,83,25;]/255;
    stimcolors_light = [240,182,211; 254 209 170;178,223,188; 255,233,166]/255;
    fullcolors = zeros(4,3,4);
    for istim = 1:4

        %tmpmap = customcolormap([0 1], [stimcolors_light(istim,:); stimcolors_dark(istim,:)],100);
        %fullcolors(:,:,istim) = tmpmap([1 35 75 100],:);

    end

    %%
    titles = {'Excitatory','Inhibitory','Uncertain'};
    for ineuron = 1:3

        groupedBarresponse(squeeze(fullcorr(:,ineuron,:)),fullcolors)
        xticklabels(neuropildata(1).dataorder);
        %%prettyaxes('ylabel','neuropil correlation','ylim',[0 1])
        title(titles{ineuron})
        legend off


    end

    istim = 3;
    ineuron = 2;
    distranges = [0 100; 100 200; 200 800];

    mydata = cell(4,3);

    for iclass = 1:4

        for idist = 1:3

         

        end

    end


    % groupedBarresponse(fullcorr)
    % xticklabels({'Excitatory','Inhibitory','Uncertain'});
    % %prettyaxes('ylabel','neuropil correlation','ylim',[0 1])
    % title(neuropildata(animalnum).dataorder{istim})
    % legend off
    % 

return

function neuropil_v_responseclass(neuropildata,istim)
%%
    fullcorr = cell(4,3);
     fulldist = cell(4,3);
    [my_annulus,responselabels,my_dist] = geteilabels(neuropildata);
    for animalnum = 1:length(neuropildata)
        dilation = 1;
        %tmpactivationmask = neuropildata(animalnum).responseclass.(responseclass)(:,istim);
        tmpemask = neuropildata(animalnum).emask;
        tmpimask = neuropildata(animalnum).imask;
        tmpmmask = neuropildata(animalnum).mmask;

        tmpdata = neuropildata(animalnum).annuluscorr(:,dilation,istim);

        % rcmask = neuropildata(animalnum).responseclass.ra(:,istim);
        % fullcorr{1,1} = cat(1,fullcorr{1,1},tmpdata(rcmask & tmpemask));
        % fullcorr{1,2} = cat(1,fullcorr{1,2},tmpdata(rcmask & tmpimask));
        % fullcorr{1,3} = cat(1,fullcorr{1,3},tmpdata(rcmask & tmpmmask));
        % 
        % rcmask = neuropildata(animalnum).responseclass.sa(:,istim);
        % fullcorr{2,1} = cat(1,fullcorr{2,1},tmpdata(rcmask & tmpemask));
        % fullcorr{2,2} = cat(1,fullcorr{2,2},tmpdata(rcmask & tmpimask));
        % fullcorr{2,3} = cat(1,fullcorr{2,3},tmpdata(rcmask & tmpmmask));
        % 
        % rcmask = neuropildata(animalnum).responseclass.na(:,istim);
        % fullcorr{3,1} = cat(1,fullcorr{3,1},tmpdata(rcmask & tmpemask));
        % fullcorr{3,2} = cat(1,fullcorr{3,2},tmpdata(rcmask & tmpimask));
        % fullcorr{3,3} = cat(1,fullcorr{3,3},tmpdata(rcmask & tmpmmask));
        % 
        % rcmask = neuropildata(animalnum).responseclass.fac(:,istim);
        % fullcorr{4,1} = cat(1,fullcorr{4,1},tmpdata(rcmask & tmpemask));
        % fullcorr{4,2} = cat(1,fullcorr{4,2},tmpdata(rcmask & tmpimask));
        % fullcorr{4,3} = cat(1,fullcorr{4,3},tmpdata(rcmask & tmpmmask));

        for iclass = 1:4

            for ineuron = 1:3
                tmpdata = my_annulus{istim,ineuron};
                tmpdist = my_dist{istim,ineuron};
                fullcorr{iclass,ineuron} = tmpdata(responselabels{istim,ineuron} == iclass);
                fulldist{iclass,ineuron} = tmpdist(responselabels{istim,ineuron} == iclass);

            end

        end

    end

    groupedBarresponse(fullcorr(:,1:2))
    xticklabels({'Excitatory','Inhibitory','Uncertain'});
    %%prettyaxes('ylabel','neuropil correlation','ylim',[0 1])
    title(neuropildata(animalnum).dataorder{istim})
    legend off

    ineuron = 2;
    distranges = [0 100; 100 200; 200 800];

    mydata = cell(4,3);

    for iclass = 1:4

        for idist = 1:3

            tmpidx = fulldist{iclass,ineuron}>distranges(idist,1) & fulldist{iclass,ineuron}<distranges(idist,2);

            if isempty(tmpidx)

                mydata{iclass,idist} = nan;
            else
                mydata{iclass,idist} = fullcorr{iclass,ineuron}(tmpidx);

            end

        end

    end

    x = cellfun(@mean,mydata);
    y = cellfun(@(x) std(x)./sqrt(length(x)),mydata);
figure
    for iclass = 1:4

        errorbar([1 2 3],x(iclass,:),y(iclass,:))
        hold on

    end

return

function plotcorr_v_annuluswidth(neuropildata,istim)
%%
    % get all neurons for each annulus width
    [annuluscorr_e,annuluscorr_i,annuluscorr_m] = deal([]);

    for animalnum = 1:length(neuropildata)
        
        % get the data
        tmpdata = neuropildata(animalnum).annuluscorr(:,:,istim);

        % get the activation and neuron subtype masks
        tmpactivationmask = repmat(neuropildata(animalnum).activationmask(:,istim),1,size(tmpdata,2));
        tmpemask = repmat(neuropildata(animalnum).emask,1,size(tmpdata,2)) & tmpactivationmask; 
        tmpimask = repmat(neuropildata(animalnum).imask,1,size(tmpdata,2)) & tmpactivationmask;
        tmpmmask = repmat(neuropildata(animalnum).mmask,1,size(tmpdata,2)) & tmpactivationmask;
        
        tmp_edata = reshape(tmpdata(tmpemask),[],size(tmpdata,2));
        tmp_idata = reshape(tmpdata(tmpimask),[],size(tmpdata,2));
        tmp_mdata = reshape(tmpdata(tmpmmask),[],size(tmpdata,2));


        annuluscorr_e = cat(1,annuluscorr_e,tmp_edata);
        annuluscorr_i = cat(1,annuluscorr_i,tmp_idata);
        annuluscorr_m = cat(1,annuluscorr_m,tmp_mdata);

    end
    

    figure,
    dilations = neuropildata(1).annuluswidth;
    for idilation = 1:size(annuluscorr_e,2)

        scatter(dilations(idilation)*ones(size(annuluscorr_e,1)),annuluscorr_e(:,idilation))
        hold on

    end

    for ineuron = 1:size(annuluscorr_e,1)

        plot(dilations,annuluscorr_e(ineuron,:))

    end

    figure,
    dilations = neuropildata(1).annuluswidth;
    for idilation = 1:size(annuluscorr_i,2)

        scatter(dilations(idilation)*ones(size(annuluscorr_i,1)),annuluscorr_i(:,idilation))
        hold on

    end
    for ineuron = 1:size(annuluscorr_i,1)

        plot(dilations,annuluscorr_i(ineuron,:))

    end



    figure,
    dilations = neuropildata(1).annuluswidth;
    for idilation = 1:size(annuluscorr_m,2)

        scatter(dilations(idilation)*ones(size(annuluscorr_m,1)),annuluscorr_m(:,idilation))
        hold on

    end
    for ineuron = 1:size(annuluscorr_m,1)

        plot(dilations,annuluscorr_m(ineuron,:))

    end




return

function plotcoor_v_dist(neuropildata,istim,dilation)
%%
    % get all neurons for each annulus width
    [annuluscorr_e,annuluscorr_i,annuluscorr_m] = deal([]);
    colors = [90,156,254;231 118 129; 157,126,213];
    colors = colors/255;
    for animalnum = 1:length(neuropildata)
        
        % get the data
        tmpdata = neuropildata(animalnum).annuluscorr(:,dilation,istim);
        tmpdist = neuropildata(animalnum).roidist;

        % get the activation and neuron subtype masks
        tmpactivationmask = neuropildata(animalnum).activationmask(:,istim);
        tmpemask = neuropildata(animalnum).emask & tmpactivationmask; 
        tmpimask = neuropildata(animalnum).imask & tmpactivationmask;
        tmpmmask = neuropildata(animalnum).mmask & tmpactivationmask;
        
        tmp_edata = [tmpdata(tmpemask), tmpdist(tmpemask)];
        tmp_idata = [tmpdata(tmpimask), tmpdist(tmpimask)];
        tmp_mdata = [tmpdata(tmpmmask), tmpdist(tmpmmask)];


        annuluscorr_e = cat(1,annuluscorr_e,tmp_edata);
        annuluscorr_i = cat(1,annuluscorr_i,tmp_idata);
        annuluscorr_m = cat(1,annuluscorr_m,tmp_mdata);

    end 


    tmpdistdata = cell(3,3);
    distanceranges = [0 100;100 200; 200 900];
    for idistrange = 1:3

        tmpidx = annuluscorr_e(:,2)>distanceranges(idistrange,1) & annuluscorr_e(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,1} = annuluscorr_e(tmpidx,1);

        tmpidx = annuluscorr_i(:,2)>distanceranges(idistrange,1) & annuluscorr_i(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,2} = annuluscorr_i(tmpidx,1);

        tmpidx = annuluscorr_m(:,2)>distanceranges(idistrange,1) & annuluscorr_m(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,3} = annuluscorr_m(tmpidx,1);

    end


    groupedBarPlotWithError(tmpdistdata(:,1:2)',colors)
    xticklabels({'<100 µm','100 - 200 µm','>200 µm'});
    %%prettyaxes('ylabel','neuropil correlation','ylim',[0 1])
    title(neuropildata(1).dataorder{istim})
    legend off
    %pbaspect([1 1.2 1])




    % organize data for stats/
    celltypes = {'E','I','M'};
    disttypes = [1 2 3];
    y = [];
    distfactor = [];
    cellfactor = {};
    for icelltype = 1:2

        for idistrange = 1:3

            tmpy = tmpdistdata{idistrange,icelltype};
            y = cat(1,y,tmpy);
            distfactor = cat(1, distfactor,repmat(disttypes(idistrange), length(tmpy),1));
            cellfactor = cat(1,cellfactor,repmat(celltypes(icelltype),length(tmpy),1));


        end

    end
    [p,tbl,stats] = anovan(y,{cellfactor,distfactor},'varnames',{'Celltype','Distance'});

    fprintf('%s Stats - Celltype: %.3e, Distance: %.3e\n',neuropildata(1).dataorder{istim},p(1),p(2) )
    %groupedviolinplot(tmpdistdata',colors)

    % group by distance bins 0-100 100-200 200+
%%
    stimname = neuropildata(animalnum).dataorder{istim};
    figure,
    scatter(annuluscorr_e(:,2),annuluscorr_e(:,1),100,colors(1,:),'filled','MarkerEdgeColor',colors(1,:) - 90/255);
    %%prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    title([stimname ' Excitatory'])
    legend off

    hold on
    scatter(annuluscorr_i(:,2),annuluscorr_i(:,1),100,colors(2,:),'filled','MarkerEdgeColor',colors(2,:) - 118/255)
    %%prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    title([stimname ' Inhibitory'])
    legend off
    %pbaspect([3 1 1])

    emdl = fitlm(annuluscorr_e(:,2),annuluscorr_e(:,1),'linear');
    y = emdl.Coefficients.Estimate(1) + emdl.Coefficients.Estimate(2)* annuluscorr_e(:,2);
    plot(annuluscorr_e(:,2), y,'LineWidth',3,'Color',colors(1,:) - 90/255)

    imdl = fitlm(annuluscorr_i(:,2),annuluscorr_i(:,1),'linear');
    y = imdl.Coefficients.Estimate(1) + imdl.Coefficients.Estimate(2)* annuluscorr_e(:,2);
    plot(annuluscorr_e(:,2), y,'LineWidth',3,'Color',colors(2,:) - 90/255)
    title(sprintf('%s Exc-p%.2e: Inh-p%.2e', stimname,emdl.Coefficients.pValue(2),imdl.Coefficients.pValue(2)))

    % figure
    % 
    % scatter(annuluscorr_m(:,2),annuluscorr_m(:,1))
    % %prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    % title([stimname ' Mixed'])
    % legend off
return


function comparisons = plotcoor_v_dist_allstim(neuropildata,dilation)
%%
    % get all neurons for each annulus width
    [annuluscorr_e,annuluscorr_i,annuluscorr_m] = deal([]);
    colors = [90,156,254;231 118 129; 157,126,213];
    colors = colors/255;
    for animalnum = 1:length(neuropildata)
        
        for istim = 1:4
            % get the data
            tmpdata = neuropildata(animalnum).annuluscorr(:,dilation,istim);
            tmpdist = neuropildata(animalnum).roidist;
    
            % get the activation and neuron subtype masks
            tmpactivationmask = neuropildata(animalnum).activationmask(:,istim);
            tmpemask = neuropildata(animalnum).emask & tmpactivationmask; 
            tmpimask = neuropildata(animalnum).imask & tmpactivationmask;
            tmpmmask = neuropildata(animalnum).mmask & tmpactivationmask;
            
            tmp_edata = [tmpdata(tmpemask), tmpdist(tmpemask)];
            tmp_idata = [tmpdata(tmpimask), tmpdist(tmpimask)];
            tmp_mdata = [tmpdata(tmpmmask), tmpdist(tmpmmask)];
    
    
            annuluscorr_e = cat(1,annuluscorr_e,tmp_edata);
            annuluscorr_i = cat(1,annuluscorr_i,tmp_idata);
            annuluscorr_m = cat(1,annuluscorr_m,tmp_mdata);

        end

    end 


    tmpdistdata = cell(3,3);
    distanceranges = [0 100;100 200; 200 900];
    for idistrange = 1:3

        tmpidx = annuluscorr_e(:,2)>distanceranges(idistrange,1) & annuluscorr_e(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,1} = annuluscorr_e(tmpidx,1);

        tmpidx = annuluscorr_i(:,2)>distanceranges(idistrange,1) & annuluscorr_i(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,2} = annuluscorr_i(tmpidx,1);

        tmpidx = annuluscorr_m(:,2)>distanceranges(idistrange,1) & annuluscorr_m(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,3} = annuluscorr_m(tmpidx,1);

    end


    groupedBarPlotWithError(tmpdistdata(:,1:2)',colors)
    xticklabels({'<100 µm','100 - 200 µm','>200 µm'});
    %%prettyaxes('ylabel','neuropil correlation','ylim',[0 1])
    title('Neuropil corr all neurons')
    legend off
    %pbaspect([1 1.2 1])

colors = [90,156,254;231 118 129; 157,126,213];
        colors = colors/255;

    figure,
    bar([1],mean(tmpdistdata{3,1}),'FaceColor',colors(1,:))
    hold on
    bar(2, mean(tmpdistdata{3,2}),'FaceColor',colors(2,:))

    errorbar([1,2],...
        [mean(tmpdistdata{3,1}),mean(tmpdistdata{3,2})],...
        [std(tmpdistdata{3,1}),std(tmpdistdata{3,2})]./sqrt([length(tmpdistdata{3,1}),length(tmpdistdata{3,2})]),...
        'linestyle','none',...
        'linewidth',1.5,...
        'color','k',...
        'capsize',10)
    xticks([1 2])
    xticklabels({'Excitatory','Inhibitory'})
    %%prettyaxes('ylabel','Correlation','ylim',[0 1])
    %pbaspect([1 2 1])
    legend off
    title('All stim conditions >200 µm')

hold off
figure

    % organize data for stats/
    celltypes = {'E','I','M'};
    disttypes = [1 2 3];
    y = [];
    distfactor = [];
    cellfactor = {};
    for icelltype = 1:2

        for idistrange = 1:3

            tmpy = tmpdistdata{idistrange,icelltype};
            y = cat(1,y,tmpy);
            distfactor = cat(1, distfactor,repmat(disttypes(idistrange), length(tmpy),1));
            cellfactor = cat(1,cellfactor,repmat(celltypes(icelltype),length(tmpy),1));


        end

    end
    [p,tbl,stats] = anovan(y,{cellfactor,distfactor},'varnames',{'Celltype','Distance'});

    [c,m,h,gnames] = multcompare(stats,'Dimension',[1,2]);

    comparisons.c = c;
    comparisons.m = m;
    comparisons.h = h;
    comparisons.gnames = gnames;
    fprintf('%s Stats - Celltype: %.3e, Distance: %.3e\n',neuropildata(1).dataorder{istim},p(1),p(2) )
    %groupedviolinplot(tmpdistdata',colors)

    % group by distance bins 0-100 100-200 200+


%%
    %stimname = neuropildata(animalnum).dataorder{istim};
    figure,
    scatter(annuluscorr_e(:,2),annuluscorr_e(:,1),100,colors(1,:),'filled','MarkerEdgeColor',colors(1,:) - 90/255);
    %prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    %title([stimname ' Excitatory'])
    legend off

    hold on
    scatter(annuluscorr_i(:,2),annuluscorr_i(:,1),100,colors(2,:),'filled','MarkerEdgeColor',colors(2,:) - 118/255)
    %prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    %title([stimname ' Inhibitory'])
    legend off
    %pbaspect([3 1 1])

    emdl = fitlm(annuluscorr_e(:,2),annuluscorr_e(:,1),'linear');
    y = emdl.Coefficients.Estimate(1) + emdl.Coefficients.Estimate(2)* annuluscorr_e(:,2);
    plot(annuluscorr_e(:,2), y,'LineWidth',3,'Color',colors(1,:) - 90/255)

    imdl = fitlm(annuluscorr_i(:,2),annuluscorr_i(:,1),'linear');
    y = imdl.Coefficients.Estimate(1) + imdl.Coefficients.Estimate(2)* annuluscorr_e(:,2);
    plot(annuluscorr_e(:,2), y,'LineWidth',3,'Color',colors(2,:) - 90/255)
    title(sprintf('All patterns Exc-p%.2e: Inh-p%.2e',emdl.Coefficients.pValue(2),imdl.Coefficients.pValue(2)))

    % figure
    % 
    % scatter(annuluscorr_m(:,2),annuluscorr_m(:,1))
    % %prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    % title([stimname ' Mixed'])
    % legend off
return

function comparisons = plotbasecoor_v_dist_allstim(neuropildata,dilation)
%%
    % get all neurons for each annulus width
    [annuluscorr_e,annuluscorr_i,annuluscorr_m] = deal([]);
    colors = [90,156,254;231 118 129; 157,126,213];
    colors = colors/255;
    for animalnum = 1:length(neuropildata)
        
        for istim = 1:4
            % get the data
            tmpdata = neuropildata(animalnum).annuluscorrbase(:,dilation,istim);
            tmpdist = neuropildata(animalnum).roidist;
    
            % get the activation and neuron subtype masks
            tmpactivationmask = neuropildata(animalnum).activationmask(:,istim);
            tmpemask = neuropildata(animalnum).emask & tmpactivationmask; 
            tmpimask = neuropildata(animalnum).imask & tmpactivationmask;
            tmpmmask = neuropildata(animalnum).mmask & tmpactivationmask;
            
            tmp_edata = [tmpdata(tmpemask), tmpdist(tmpemask)];
            tmp_idata = [tmpdata(tmpimask), tmpdist(tmpimask)];
            tmp_mdata = [tmpdata(tmpmmask), tmpdist(tmpmmask)];
    
    
            annuluscorr_e = cat(1,annuluscorr_e,tmp_edata);
            annuluscorr_i = cat(1,annuluscorr_i,tmp_idata);
            annuluscorr_m = cat(1,annuluscorr_m,tmp_mdata);

        end

    end 


    tmpdistdata = cell(3,3);
    distanceranges = [0 100;100 200; 200 900];
    for idistrange = 1:3

        tmpidx = annuluscorr_e(:,2)>distanceranges(idistrange,1) & annuluscorr_e(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,1} = annuluscorr_e(tmpidx,1);

        tmpidx = annuluscorr_i(:,2)>distanceranges(idistrange,1) & annuluscorr_i(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,2} = annuluscorr_i(tmpidx,1);

        tmpidx = annuluscorr_m(:,2)>distanceranges(idistrange,1) & annuluscorr_m(:,2)<=distanceranges(idistrange,2);
        tmpdistdata{idistrange,3} = annuluscorr_m(tmpidx,1);

    end


    groupedBarPlotWithError(tmpdistdata(:,1:2)',colors)
    xticklabels({'<100 µm','100 - 200 µm','>200 µm'});
    %prettyaxes('ylabel','neuropil correlation','ylim',[0 1])
    title('Neuropil corr all neurons')
    legend off
    %pbaspect([1 1.2 1])

colors = [90,156,254;231 118 129; 157,126,213];
        colors = colors/255;

    figure,
    bar([1],mean(tmpdistdata{3,1}),'FaceColor',colors(1,:))
    hold on
    bar(2, mean(tmpdistdata{3,2}),'FaceColor',colors(2,:))

    errorbar([1,2],...
        [mean(tmpdistdata{3,1}),mean(tmpdistdata{3,2})],...
        [std(tmpdistdata{3,1}),std(tmpdistdata{3,2})]./sqrt([length(tmpdistdata{3,1}),length(tmpdistdata{3,2})]),...
        'linestyle','none',...
        'linewidth',1.5,...
        'color','k',...
        'capsize',10)
    xticks([1 2])
    xticklabels({'Excitatory','Inhibitory'})
    %prettyaxes('ylabel','Correlation','ylim',[0 1])
    %pbaspect([1 2 1])
    legend off
    title('All stim conditions >200 µm')

hold off
figure

    % organize data for stats/
    celltypes = {'E','I','M'};
    disttypes = [1 2 3];
    y = [];
    distfactor = [];
    cellfactor = {};
    for icelltype = 1:2

        for idistrange = 1:3

            tmpy = tmpdistdata{idistrange,icelltype};
            y = cat(1,y,tmpy);
            distfactor = cat(1, distfactor,repmat(disttypes(idistrange), length(tmpy),1));
            cellfactor = cat(1,cellfactor,repmat(celltypes(icelltype),length(tmpy),1));


        end

    end
    [p,tbl,stats] = anovan(y,{cellfactor,distfactor},'varnames',{'Celltype','Distance'});

    [c,m,h,gnames] = multcompare(stats,'Dimension',[1,2]);

    comparisons.c = c;
    comparisons.m = m;
    comparisons.h = h;
    comparisons.gnames = gnames;
    fprintf('%s Stats - Celltype: %.3e, Distance: %.3e\n',neuropildata(1).dataorder{istim},p(1),p(2) )
    %groupedviolinplot(tmpdistdata',colors)

    % group by distance bins 0-100 100-200 200+


%%
    %stimname = neuropildata(animalnum).dataorder{istim};
    figure,
    scatter(annuluscorr_e(:,2),annuluscorr_e(:,1),100,colors(1,:),'filled','MarkerEdgeColor',colors(1,:) - 90/255);
    %prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    %title([stimname ' Excitatory'])
    legend off

    hold on
    scatter(annuluscorr_i(:,2),annuluscorr_i(:,1),100,colors(2,:),'filled','MarkerEdgeColor',colors(2,:) - 118/255)
    %prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    %title([stimname ' Inhibitory'])
    legend off
    %pbaspect([3 1 1])

    emdl = fitlm(annuluscorr_e(:,2),annuluscorr_e(:,1),'linear');
    y = emdl.Coefficients.Estimate(1) + emdl.Coefficients.Estimate(2)* annuluscorr_e(:,2);
    plot(annuluscorr_e(:,2), y,'LineWidth',3,'Color',colors(1,:) - 90/255)

    imdl = fitlm(annuluscorr_i(:,2),annuluscorr_i(:,1),'linear');
    y = imdl.Coefficients.Estimate(1) + imdl.Coefficients.Estimate(2)* annuluscorr_e(:,2);
    plot(annuluscorr_e(:,2), y,'LineWidth',3,'Color',colors(2,:) - 90/255)
    title(sprintf('All patterns Exc-p%.2e: Inh-p%.2e',emdl.Coefficients.pValue(2),imdl.Coefficients.pValue(2)))

    % figure
    % 
    % scatter(annuluscorr_m(:,2),annuluscorr_m(:,1))
    % %prettyaxes('ylabel','Neuropil correlation','xlabel','Distance, µm')
    % title([stimname ' Mixed'])
    % legend off
return

function groupedviolinplot(data,colors)

    numGroups = size(data, 1);
    numConditions = size(data, 2);
    
    % Define colors for each group
    %groupColors = colors;
    
    % Bar width and spacing
    barWidth = 0.5;  % Adjust as needed
    barSpacing = 0.15;  % Adjust as needed
    
    % X positions for the groups
    xPositions = 1:2:numConditions*2;
    
    % Create grouped bar plot with error bars
    figure;
    hold on;
    
    for i = 1:numGroups
        %for j = 1:numConditions
            % Calculate x position for each bar in the group
            xPos = [xPositions(i) + (i - 1) * (barWidth + barSpacing), xPositions(i) + (i ) * (barWidth + barSpacing),xPositions(i) + (i + 1) * (barWidth + barSpacing)]
            
            % Calculate mean and standard error for the current cell
            % meanValue = mean(data{i, j});
            % stdError = std(data{i, j}) / sqrt(length(data{i, j}));
            
            % Plot the bar with error bar
            %bar(xPos, meanValue, 'BarWidth', barWidth, 'FaceColor', groupColors(i, :));
            %errorbar(xPos, meanValue, stdError, 'k', 'LineWidth', 1.5);  % 'k' for black color

            violin(data(i,:),'x',xPos,'facecolor',colors,'facealpha',0.7,'bw',0.05);
            for j = 1:3

                scatter(xPos(j)*ones(length(data{i,j})),data{i,j})

            end

        %end
    end
    
    % Customize the plot as needed
    xticks([1.5:2:7.5])
    xlabel('Conditions');
    ylabel('Mean Values');
    title('Grouped Bar Plot with Error Bars');
    legend({'Group 1', 'Group 2', 'Group 3'});  % Add group labels
    axis tight
    
    hold off;


return

function meancorrelation(neuropildata)
%%
    fullcorr = cell(3,4);

    for animalnum = 1:length(neuropildata)
        dilation = 1;
        tmpactivationmask = neuropildata(animalnum).activationmask;
        tmpemask = neuropildata(animalnum).emask;
        tmpimask = neuropildata(animalnum).imask;
        tmpmmask = neuropildata(animalnum).mmask;


        for istim = 1:4
            
            tmpdata = neuropildata(animalnum).annuluscorr(:,dilation,istim);
            
            fullcorr{1,istim} = cat(1,fullcorr{1,istim}, tmpdata(tmpactivationmask(:,istim) & tmpemask));
            fullcorr{2,istim} = cat(1,fullcorr{2,istim}, tmpdata(tmpactivationmask(:,istim)  & tmpimask));
            fullcorr{3,istim} = cat(1,fullcorr{3,istim}, tmpdata(tmpactivationmask(:,istim)  & tmpmmask));


        end


    end

    groupedBarPlotWithError(fullcorr(1:2,:))
    xticklabels(neuropildata(1).dataorder);
    %prettyaxes('ylabel','neuropil correlation','ylim',[0 1])
    legend off

    % prepare the stats
    stimnames = neuropildata(1).dataorder;
    cells = {'E','I','M'};
    y = [];
    cellfactor = {};
    stimfactor = {};
    for istim = 1:4

        for icell = 1:2

            tmpy = fullcorr{icell,istim};
            y = cat(1,y,tmpy);

            cellfactor = cat(1,cellfactor,repmat(cells(icell),length(tmpy),1));
            stimfactor = cat(1,stimfactor,repmat(stimnames(istim),length(tmpy),1));

        end

    end

    [p,t,stats] = anovan(y,{cellfactor,stimfactor},'varnames',{'Cell','Stim'});
    fprintf('Stats: Cellfactor: %.3e, Stimfactor: %.3e', p(1),p(2));
    [c,m,h,gnames] = multcompare(stats,'Dimension',[1,2]);

    pvals = zeros(1,4);
    for istim = 1:4

        [h,p] = ttest2(fullcorr{1,istim},fullcorr{2,istim});
        pvals(istim) = p;

    end
    %pvals = bonf_holm(pvals,0.05);

return

function groupedBarPlotWithError(data,colors)
    % Input:
    %   data: 3x4 cell array representing 3 groups for 4 different conditions
    %         Each cell contains a different number of data points
    

    if ~exist('colors','var')
    
        colors = [90,156,254;231 118 129; 157,126,213];
        colors = colors/255;

    end

    % Number of groups and conditions
    numGroups = size(data, 1);
    numConditions = size(data, 2);
    
    % Define colors for each group
    groupColors = colors;
    
    % Bar width and spacing
    barWidth = 0.3;  % Adjust as needed
    barSpacing = 0.08;  % Adjust as needed
    
    % X positions for the groups
    xPositions = 1:numConditions;
    
    % Create grouped bar plot with error bars
    figure;
    hold on;
    
    for i = 1:numGroups
        for j = 1:numConditions
            % Calculate x position for each bar in the group
            xPos = xPositions(j) + (i - 1) * (barWidth + barSpacing);
            
            % Calculate mean and standard error for the current cell
            meanValue = mean(data{i, j});
            stdError = std(data{i, j}) / sqrt(length(data{i, j}));
            
            % Plot the bar with error bar
            bar(xPos, meanValue, 'BarWidth', barWidth, 'FaceColor', groupColors(i, :));
            errorbar(xPos, meanValue, stdError, 'k', 'LineWidth', 1.5);  % 'k' for black color
        end
    end
    
    % Customize the plot as needed
    xticks([1.25:4.25])
    xlabel('Conditions');
    ylabel('Mean Values');
    title('Grouped Bar Plot with Error Bars');
    legend({'Group 1', 'Group 2', 'Group 3'});  % Add group labels
    
    hold off;
return

function groupedBarresponse(data,colors)
    % Input:
    %   data: 3x4 cell array representing 3 groups for 4 different conditions
    %         Each cell contains a different number of data points
    
%%
    if ~exist('colors','var')
    
        colors = [90,156,254;231 118 129; 157,126,213];
        colors = colors/255;

        fullcolors = zeros(4,3,3);
        for icell = 1:3
        
            %tmpmap = customcolormap([0  1],[colors(icell,:) - min(colors(icell,:));colors(icell,:) ],100);
            %fullcolors(:,:,icell) = tmpmap([1,30,60,100],:);

        end

    else

        fullcolors = colors

    end
%%
    % Number of groups and conditions
    numGroups = size(data, 1);
    numConditions = size(data, 2);
    
    % Define colors for each group
   
    
    % Bar width and spacing
    barWidth = 0.3;  % Adjust as needed
    barSpacing = 0.15;  % Adjust as needed
    
    % X positions for the groups
    xPositions = 1:2:numConditions*2;
    
    % Create grouped bar plot with error bars
    figure;
    hold on;
    
    for i = 1:numGroups
        for j = 1:numConditions
            % Calculate x position for each bar in the group
            xPos = xPositions(j) + (i - 1) * (barWidth + barSpacing);
            
            % Calculate mean and standard error for the current cell
            meanValue = mean(data{i, j});
            stdError = std(data{i, j}) / sqrt(length(data{i, j}));
            
            % Plot the bar with error bar
            bar(xPos, meanValue, 'BarWidth', barWidth, 'FaceColor', fullcolors(i, :,j));
            errorbar(xPos, meanValue, stdError, 'k', 'LineWidth', 1.5);  % 'k' for black color
        end
    end
    
    % Customize the plot as needed
    xticks([1.5:2:7.5])
    xlabel('Conditions');
    ylabel('Mean Values');
    title('Grouped Bar Plot with Error Bars');
    legend({'Group 1', 'Group 2', 'Group 3'});  % Add group labels
    
    hold off;
return


function [my_annulus,responselabels,my_dist] = geteilabels(neuropildata)
%%
    frate = 0.034; %frame rate
    x_label = [1:1800]*frate; %time axis
    onset_start = 294; %approximate frame that ICMS started
    onset_stop = floor(onset_start+3/frate);
    offset_start = 1177;  %approximate frame ICMS ended
    offset_stop = floor(offset_start+3/frate);
    num_animals = 4; %total number of animals
    animals = {'EIM08', 'EIM09', 'VGAT2','VGAT5'}; %animal IDs
    train_types = {'10 Hz', '10 Hz Burst', '100 Hz', 'TBS'}; %different ICMS trains used 

    % organize data like chris did
    my_data = cell(4,3);
    my_annulus = cell(4,3);
    my_dist = cell(4,3);

    for animalnum = 1:4

        tmpemask = neuropildata(animalnum).emask;
        tmpimask = neuropildata(animalnum).imask;
        tmpmmask = neuropildata(animalnum).mmask;
        for istim = 1:4
            tmpactivationmask = neuropildata(animalnum).activationmask(:,istim);
    
            tmpdata = neuropildata(animalnum).fijicalcium(:,:,istim);
            tmpannulus = neuropildata(animalnum).annuluscorr(:,1,istim);
            tmpannulus(~tmpactivationmask) = nan;
            tmpdata(:,~tmpactivationmask) = nan;
            tmpdist = neuropildata(animalnum).roidist;
    
            my_data{istim,1} = cat(1,my_data{istim,1},tmpdata(:,tmpemask)');
            my_data{istim,2} = cat(1,my_data{istim,2},tmpdata(:,tmpimask)');
            my_data{istim,3} = cat(1,my_data{istim,3},tmpdata(:,tmpmmask)');

            my_annulus{istim,1} = cat(1,my_annulus{istim,1},tmpannulus(tmpemask));
            my_annulus{istim,2} = cat(1,my_annulus{istim,2},tmpannulus(tmpimask));
            my_annulus{istim,3} = cat(1,my_annulus{istim,3},tmpannulus(tmpmmask));

            my_dist{istim,1} = cat(1,my_dist{istim,1},tmpdist(tmpemask));
            my_dist{istim,2} = cat(1,my_dist{istim,2},tmpdist(tmpimask));
            my_dist{istim,3} = cat(1,my_dist{istim,3},tmpdist(tmpmmask));

        end

    end

    my_error = cellfun(@(x) std(x(:,1:onset_start),[],2),my_data,'UniformOutput',false);



    labels = labels_during_stim(my_data,my_error);
    
    numneurons = cellfun(@length,my_annulus);

    responselabels = cell(size(my_annulus));

    for istim = 1:4
    
        for iclass = 1:3

            responselabels{istim,iclass} = labels{istim}(iclass,1:numneurons(1,iclass))';

        end

    end



    
    

return
