% nsdStim_expand.m
%
% associated with the following publication: Roth, ZN, and Merriam, EP (2023).
% Representations in human primary visual cortex drift over time
% DOI:
%
%   usage: nsdStim_expand()
%   by: zvi roth
%   date: 3/10/2022
%   purpose: run natural scene stimuli through the models, get filter energy responses,
%   and save output
%   creates files used by: prfSampleModel_expand.m

% uses the steerable pyramid: https://github.com/elimerriam/stimulusVignetting


%expanded model for Representational Drift. Includes 2 additional filters,
%for highest and lowest spatial frequencies, but with no orientation
%tuning.
%also uses sqrt instead of filter energy.

close all
clear all

interpImgSize = 714;
backgroundSize = 1024;
imgScaling = 0.5;

normResp=0;
pyramidfolder = '/misc/data18/rothzn/nsd/stimuli/pyramid_expand/';%to save model outputs

%%

% construct quad frequency filters
numOrientations = 8;
bandwidth = 1;
dims = [backgroundSize backgroundSize];
dims = dims*imgScaling;
numLevels = maxLevel(dims,bandwidth);
[freqResps, pind] = makeSteerFRs(dims, numLevels, numOrientations, bandwidth);
%%
stimfolder = '/misc/data18/rothzn/nsd/stimuli/';
stimfilename = fullfile(stimfolder,'nsd_stimuli.hdf5');%[3 1360 714 220]
stiminfo = h5info(stimfilename);

imgSizeX = stiminfo.Datasets.Dataspace.Size(2);%1360;
imgSizeY = stiminfo.Datasets.Dataspace.Size(3);%714;
numImgs = stiminfo.Datasets.Dataspace.Size(4);%220

nsdfolder = '/misc/data18/rothzn/nsd/';
nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
allImgs = nsdDesign.sharedix; %indices of the shared 1000 images


for isub=1:8
    
    
    allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
    allImgs = unique(allImgs);
    %%
    backgroundColor(1,1,:) = uint8([127,127,127]);
    
    
    nImgs = 1;

    fixPoint(1,1,:) = [255 0 0];
    
    iimg=0;
    for imgNum=allImgs
        iimg = iimg+1
        
        pyramidfilename = ['pyrImg' num2str(imgNum) '.mat'];
        if ~isfile(fullfile(pyramidfolder, pyramidfilename))%if file exists already no need to remake it
            origImg = h5read(stimfilename,'/imgBrick/',[1 1 1 imgNum],[3 imgSizeX imgSizeY nImgs]);
            origImg = double(origImg);
            origImg = permute(origImg,[3 2 1]);%[425,425,3]
            
            [Xq, Yq] = meshgrid(linspace(1,imgSizeX, interpImgSize), linspace(1,imgSizeY, interpImgSize));
            for irgb=1:3
                interpImg(:,:,irgb) = interp2(squeeze(origImg(:,:,irgb)), Xq, Yq);
            end
            %add red semi-transparent fixation point
            
            interpImg(interpImgSize/2-8:interpImgSize/2+8,interpImgSize/2-8:interpImgSize/2+8,:) = ...
                (interpImg(interpImgSize/2-8:interpImgSize/2+8,interpImgSize/2-8:interpImgSize/2+8,:) + ...
                repmat(fixPoint,17,17,1))/2;
            
            %%add background
            bigImg = repmat(backgroundColor,backgroundSize,backgroundSize,1);
            bigImg(1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2, 1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2,:) = interpImg(:,:,:);
            
            %change to grayscale by averaging across RGB channels
            bigImg = mean(bigImg,3);

            %DOWNSAMPLE
            bigImg = imresize(bigImg,imgScaling);
            %% pass image through steerable pyramid
            
            [pyr, pind] = buildSteerBands(bigImg, freqResps);
            sumOri = cell(numLevels+2,1);
            modelOri = cell(numLevels+2,1);
            for ilev = 1:numLevels
                % loop over levels and orientations of the pyramid
                % initialize output
                sumOri{ilev}(:,:) = zeros(dims(1), dims(2));
                modelOri{ilev} = zeros(numOrientations, dims(1), dims(2));
                for orientation = 1:numOrientations
                    if normResp
                        nEnergies = normEnergies(pyr,pind,numOrientations,0.1);
                        thisBand = abs(accessSteerBand(nEnergies,pind,numOrientations,ilev,orientation));
                    else
%                         thisBand = abs(accessSteerBand(pyr, pind, numOrientations,ilev, orientation)).^2;
                        thisBand = abs(accessSteerBand(pyr, pind, numOrientations,ilev, orientation));
                    end
                    sumOri{ilev}(:,:) = sumOri{ilev}(:,:) + thisBand;
                    modelOri{ilev}(orientation,:,:) = thisBand;
                end
            end
            % ADD LOWEST SPATIAL FREQUENCY
            thisBand = abs(pyrLow(pyr,pind));
%             thisBand = abs(pyrLow(freqResps,pind)).^2;
            sumOri{numLevels+1} = thisBand;
            modelOri{numLevels+1} = thisBand;
            
            % ADD HIGHEST SPATIAL FREQUENCY
            thisBand = abs(pyrHi(pyr,pind));
%           thisBand = abs(pyrHi(freqResps,pind)).^2;
            sumOri{numLevels+2} = thisBand;
            modelOri{numLevels+2} = thisBand;

            save(fullfile(pyramidfolder, pyramidfilename), 'interpImgSize','backgroundSize','imgScaling',...
                'numOrientations','bandwidth','dims','bigImg','sumOri','modelOri','numLevels','normResp');
        end
    end
end

