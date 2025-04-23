function plotneuropilanalysis(neuropildata)


    % get the activated cells
    % get the e/i/other mask
    % get the fascilitating depressing neurons

    % plot the neuropilcorr across the patterns for all activated neurons
    % plot the eneuropill corr across patterns for e/i neurons

    % plot the neuropil corr across patterns for facilitating depressing


    % get the activation mask for all cels all animals
    basetime = 10;
    stimtime = 30;
    threshmult = 3;
    rmsmult = 1;
    fullactivationmask = cell(4,1);
    [emask,imask,othermask] = deal(cell(4,1));
    waitfor(msgbox('Select the github base directory'))
    basedir = uigetdir();

    for animalnum = 1:length(neuropildata)
        
        % activation mask
        tmpactivationmask = false(length(neuropildata(animalnum).annuluscorr),4);
        for istim = 1:4
         tmpactivationmask(:,istim)= getactivationmask(neuropildata(animalnum).cellcalcium(:,:,istim),...
            neuropildata(animalnum).frate,...
            basetime,stimtime,threshmult,true);

        end
        fullactivationmask{animalnum} = tmpactivationmask;

        % ei mask
        [emask{animalnum},imask{animalnum},othermask{animalnum}] = geteimask(basedir,...
            neuropildata(animalnum).id,...
            rmsmult,basetime,neuropildata(animalnum).frate);

        % fascilitating/depressing mask (not yet set up)


    end


    %%

        eneurons = getcellmask(fullactivationmask,emask);
        ineurons = getcellmask(fullactivationmask,imask);

        ecorr = organizecellcorr(neuropildata,eneurons,'annuluscorr');
        icorr = organizecellcorr(neuropildata,ineurons,'annuluscorr');

    violin2groups(ecorr,icorr)

    %%

function cellmask = getcellmask(activationmask,eimask)
    
    if isempty(eimask)
        cellmask = activationmask;
    else
        cellmask = cell(size(activationmask));
        for animalnum = 1:length(activationmask)

            cellmask{animalnum} = activationmask{animalnum} & eimask{animalnum};

        end
    end

function cellcorr = organizecellcorr(neuropildata,cellmask,corrfield)


    cellcorr = cell(4,1);

    for animalnum = 1:length(cellmask)

        for istim = 1:4

            cellcorr{istim} = cat(1,cellcorr{istim},...
                neuropildata(animalnum).(corrfield)(cellmask{animalnum}(:,istim),istim));

        end

    end


function violin2groups(group1corr,group2corr)

    colors = [0.1 0.1 0.1;
        0.9 0.1 0.1];

    gap = 0.35;
    positions = [];
    
    spacing = [1 3 5 7];
    for igrouping = 1:4
        
        tmpposition = [spacing(igrouping)-gap, spacing(igrouping)+gap];
        positions = [positions,tmpposition];
        
    end
    
    fullcorr = cat(1,group1corr',group2corr');
    positions = [positions(1:2:8);positions(2:2:8)];
    figure
    ip = 1;
    for igroup = 1:2

        for istim = 1:4

            v(ip) = violinplot(fullcorr{igroup,istim},positions(igroup,istim),[],'Width',0.3,'ShowMean',true,'Boxwidth',0.03,'EdgeColor',colors(1,:),'ViolinColor',colors(1,:),'ViolinAlpha',0.6,'ShowData',true);
            ip = ip+1;
        end

    end

    for p = 5:8

        v(p).ViolinColor = colors(2,:);

    end

    for p = 1:8
        
        v(p).BoxPlot.LineWidth = 1;
        v(p).WhiskerPlot.LineWidth = 2;
        v(p).WhiskerPlot.LineStyle = '-';
        v(p).MedianPlot.SizeData = 50;
        v(p).MedianPlot.Marker = 'diamond';
    %gp_violin(p).MedianPlot.Marker = 'o';
        v(p).ViolinPlot.LineStyle = 'none';

    end



function violinonegroup(cellcorr)
    hexcolors = {'#7A306C','#8E8DBE','#A9E4EF','#81F495'};
    colors = zeros(4,3);
    for icolor = 1:4

        colors(icolor,:) = hex2rgb(hexcolors{icolor},1);

    end

    figure,
    for istim = 1:4

       v(istim) = violinplot(cellcorr{istim},istim);

    end
    
    % pretty
    for p = 1:4
        v(p).ViolinColor = colors(p,:);
        v(p).BoxPlot.LineWidth = 1;
        v(p).WhiskerPlot.LineWidth = 2;
        v(p).WhiskerPlot.LineStyle = '-';
        v(p).MedianPlot.SizeData = 50;
        v(p).MedianPlot.Marker = 'diamond';
    %gp_violin(p).MedianPlot.Marker = 'o';
        v(p).ViolinPlot.LineStyle = 'none';

    end


