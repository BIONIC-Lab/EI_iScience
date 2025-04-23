function [roiindices,annulusindices,annulusmask] = getneuropilannulus(img,ellipseparams,extensionMults,showmasks)

%% 
if ~exist("showmasks",'var')
    showmasks = false;
end

% get image size and the original ellipse params
[imgHeight,imgWidth] = size(img);
centerX = ellipseparams(:,1);
centerY = ellipseparams(:,2);
majorAxis = ellipseparams(:,3)/2;
minorAxis = ellipseparams(:,4)/2;
diameter = mean(ellipseparams(:,3:4),2)/2;
angle = ellipseparams(:,5);
numrois = size(ellipseparams,1);

% Generate the binary mask for each roi to make sure annulus doesnt overlap
blackimg = zeros(imgHeight,imgWidth);
fullmask = false(imgHeight,imgWidth);
numberedmask = nan(imgHeight,imgWidth);
annulusmask = nan(imgHeight,imgWidth,numrois,length(extensionMults));


for iroi = 1:numrois
    
    % make elliptical roi and make a mask of roi
    roi(iroi) = images.roi.Ellipse([],'Center',[centerX(iroi),centerY(iroi)],'Semiaxes',[majorAxis(iroi),minorAxis(iroi)],'RotationAngle',angle(iroi));
    tmpmask = createMask(roi(iroi),blackimg);

    % dilate the roi for the annulus in 3 different multiples excluding
    % a ring half the diameter of the neuron 
    tmpdiameter = diameter(iroi);
    tmpannulus = zeros(imgHeight,imgWidth,length(extensionMults));
    for idilation = 1:length(extensionMults)

        tmpextension = round(extensionMults(idilation)*tmpdiameter);
        tmpannulus(:,:,idilation) = double(imdilate(tmpmask,strel('disk',tmpextension)) & ~imdilate(tmpmask,strel('disk',round(tmpdiameter*0.2))));

    end

    % set the roi to the number of the roi for later
    numberedmask(tmpmask & isnan(numberedmask)) = iroi;
    
    %set the annulus
    tmpannulus(tmpannulus==1) = iroi;
    for idilation =1:length(extensionmults)

        annulusmask(:,:,iroi,:) = tmpannulus;

    end
    fullmask = fullmask | createMask(roi(iroi),blackimg);


    
    

end

annulusmask(~isnan(repmat(numberedmask,1,1,numrois,length(extensionmults)))) = nan;
annulusmask(annulusmask ==0) = nan;

% get the indices
[roiindices] = cell(numrois,1);
annulusindices = cell(numrois,length(extensionmults));
for iroi = 1:numrois

    roiindices{iroi} = find(numberedmask == iroi);
    
    for idilation = 1:length(extensionmults)

        annulusindices{iroi,idilation} = find(annulusmask(:,:,iroi,idilation)==iroi);

    end

end

% % generate the annulus masks and make sure they dont overlap with existing
% % ones
% [annulusMasks,roimasks] = deal(false(imageHeight,imageWidth,numrois));
% [neuropilindices,roiindices] = deal(cell(numrois,1));
% for i = 1:numrois
% 
%     originalEllipseMask = createMask(roi(i), img); % create the binary mask for the roi
%     roimasks(:,:,i) = originalEllipseMask; % save the roi mask
% 
%     % Dilate the original ellipse mask to create the outer ellipse mask
%     outerEllipseMask = imdilate(originalEllipseMask, strel('disk', extensionMults));
% 
%     % Check for overlap with previously generated ellipses and annulus masks
%     annulusMask = outerEllipseMask & ~fullellipsemask;
%     for j = 1:(i-1)
%         overlapEllipse = annulusMask & annulusMasks(:,:,j);
%         overlapOriginal = annulusMask & ~annulusMasks(:,:,j) & originalEllipseMask;
%         if any(overlapEllipse(:)) || any(overlapOriginal(:))
%             annulusMask(overlapEllipse | overlapOriginal) = 0;
%         end
%     end
% 
%     % Store the annulus mask
%     annulusMasks(:,:,i) = annulusMask;
% 
%     % convert the masks into indices
%     neuropilindices{i} = find(annulusMask);
%     roiindices{i} = find(originalEllipseMask);
% 
% end

% Display the annulus masks (optional)
if showmasks
    figure
    imagesc(sum(roimasks,3).*img),colormap gray, axis square
    figure
    imagesc(sum(annulusMasks,3).*img),colormap gray, axis square
end


% convert the mask information into indices to grab.
