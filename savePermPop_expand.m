% savePermPop_expand.m
%
% associated with the following publication: Roth, ZN, and Merriam, EP (2023).
% Representations in human primary visual cortex drift over time
% DOI:
%
%   usage: savePermPop_expand(1,0,1000)
%   by: zvi roth
%   date: 3/10/2022
%   purpose:  permute the session order and compute measures of population representational drift
%   uses files created by: simPopResponse_expand.m
%   creates files used by: 

function res = savePermPop_expand(rois,toZscore,nperms,r2thresh,fixedFirst)
tic
if ieNotDefined('rois'), rois = [1]; end 
if ieNotDefined('toZscore'), toZscore = 0; end%
if ieNotDefined('nperms'), nperms = 2; end
if ieNotDefined('r2thresh'), r2thresh = 0; end
if ieNotDefined('fixedFirst'), fixedFirst = 0; end

fixedFirstStr='';
if fixedFirst
    fixedFirstStr = '_fixedFirst_';
end
rng(1)

saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand','/');
if ~isfolder(saveFolder)
    saveFolder = ['/misc/data18/rothzn/nsd/repDrift_expand/'];
end


corrTypeBetweenRDMs = 'Spearman';%'Pearson'; %for comparing RDMs across sessions
distTypeRDM = 'correlation';%'correlation', 'euclidean', 'mahalanobis', 'spearman'. For creating RDM
nrois = 4;
subjects = 1:8;%[5:8];
nsubjects = length(subjects);

nimg = 100;
nsessions=30;
maxSessions = nsessions;
minSessions=nsessions;

r2threshStr = '';
if r2thresh>0
   r2threshStr = ['r2thresh' num2str(r2thresh,'%4.2f')];
end
errorbarColor = [0.2 0.2 0.2];
surfaceAlpha = 0.1;
linewidth = 2;



zscoreStr='';
if toZscore==1
    zscoreStr = '_zscore';
elseif toZscore==2
    zscoreStr = '_zeroMean';
elseif toZscore==3
    zscoreStr = '_equalStd';
elseif toZscore==4
    zscoreStr = '_zeroROImean';
end

imgCorrMat = NaN(nrois,length(subjects),nimg,nsessions,nsessions);
imgCorrMatOri = NaN(nrois,length(subjects),nimg,nsessions,nsessions);
imgCorrVec = NaN(nrois,length(subjects),nimg,nsessions*(nsessions-1)/2);
imgCorrVecOri = NaN(nrois,length(subjects),nimg,nsessions*(nsessions-1)/2);
sessCorrMat = NaN(nrois,length(subjects),nsessions,nimg,nimg);
sessCorrMatOri = NaN(nrois,length(subjects),nsessions,nimg,nimg);
sessCorrVec = NaN(nrois,length(subjects),nsessions,nimg*(nimg-1)/2);
sessCorrVecOri = NaN(nrois,length(subjects),nsessions,nimg*(nimg-1)/2);
betweenSessCorr = NaN(nrois,length(subjects),nsessions,nsessions);
betweenSessCorrOri = NaN(nrois,length(subjects),nsessions,nsessions);

for isub=1:length(subjects)
    isub
    load(fullfile(saveFolder,['regressSessCombineROI_sub' num2str(subjects(isub)) zscoreStr '.mat']),'allRoiPrf','nsplits');
    nsessions = nsplits-1; %maxnumber of  sessions. nsplits includes the mean, and is <41 for some subjects
    subSessions(isub,1:nsessions) = ones;

    for iroi=rois%1:nrois
        load(fullfile(saveFolder,['simPopResp_v' num2str(iroi) '_sub' num2str(subjects(isub)) zscoreStr '.mat']),...
            'voxSessResp','voxSessRespOri', 'simImgs','nsessions');
        subNumSessions(isub) = nsessions;

        goodVox = allRoiPrf{iroi}.r2>r2thresh;
        r2GoodVox(iroi,isub) = r2thresh;
        numGoodVox(iroi,isub) = sum(goodVox);
        
        for iimg=1:nimg
            %for each image, a correlation matrix between sessions
            imgCorrMat(iroi,isub,iimg,1:nsessions,1:nsessions) = corr(squeeze(voxSessResp(goodVox,:,iimg)));
            imgCorrMatOri(iroi,isub,iimg,1:nsessions,1:nsessions) = corr(squeeze(voxSessRespOri(goodVox,:,iimg)));
            imgCorrVec(iroi,isub,iimg,1:nsessions*(nsessions-1)/2) = pdist(squeeze(voxSessResp(goodVox,:,iimg))',distTypeRDM);
            imgCorrVecOri(iroi,isub,iimg,1:nsessions*(nsessions-1)/2) = pdist(squeeze(voxSessRespOri(goodVox,:,iimg))',distTypeRDM);
        end
        
        for isess=1:nsessions
            %RDM of all images in a single session
            sessCorrMat(iroi,isub,isess,:,:) = corr(squeeze(voxSessResp(goodVox,isess,:)));
            sessCorrMatOri(iroi,isub,isess,:,:) = corr(squeeze(voxSessRespOri(goodVox,isess,:)));
            sessCorrVec(iroi,isub,isess,:) = pdist(squeeze(voxSessResp(goodVox,isess,:))',distTypeRDM);
            sessCorrVecOri(iroi,isub,isess,:) = pdist(squeeze(voxSessRespOri(goodVox,isess,:))',distTypeRDM);
            
        end
        %correlate between RDMs of different sessions
        betweenSessCorr(iroi,isub,1:minSessions,1:minSessions) = corr(squeeze(sessCorrVec(iroi,isub,1:minSessions,1:minSessions))','Type',corrTypeBetweenRDMs);
        betweenSessCorrOri(iroi,isub,1:minSessions,1:minSessions) = corr(squeeze(sessCorrVecOri(iroi,isub,1:minSessions,1:minSessions))','Type',corrTypeBetweenRDMs);
        
    end
end
minSessions = min(subNumSessions);

%average across images
avgImgCorrMat = squeeze(mean(imgCorrMat,3));
avgImgCorrMatOri = squeeze(mean(imgCorrMatOri,3));

%initialize permutation matrices of average across images.
imgCorrMatPerm = NaN(nrois,length(subjects),nperms,minSessions,minSessions);
imgCorrMatOriPerm = NaN(nrois,length(subjects),nperms,minSessions,minSessions);

%initialize permutation matrices of correlations between RDMs.
betweenSessCorrPerm = NaN(nrois,length(subjects),nperms,minSessions,minSessions);
betweenSessCorrOriPerm = NaN(nrois,length(subjects),nperms,minSessions,minSessions);

%%
distMatrix = toeplitz(0:nsessions-1);

betweenSessImg = NaN(nrois,length(subjects),minSessions-1);
betweenSessImgOri = NaN(nrois,length(subjects),minSessions-1);
betweenSessImgPerm = NaN(nrois,length(subjects),nperms,minSessions-1);
betweenSessImgOriPerm = NaN(nrois,length(subjects),nperms,minSessions-1);
betweenSessDist = NaN(nrois,length(subjects),minSessions-1);
betweenSessDistOri = NaN(nrois,length(subjects),minSessions-1);
betweenSessDistPerm = NaN(nrois,length(subjects),nperms,minSessions-1);
betweenSessDistOriPerm = NaN(nrois,length(subjects),nperms,minSessions-1);

for iroi=rois%1:nrois
    %correlations between population responses
    for isub=1:length(subjects)
        nsessions = subNumSessions(isub);
        if isub>1%use subject 1 randomized orders
            permOrders{isub} = permOrders{1};
        else%randomize for subject 1
            for iperm=1:nperms
                if fixedFirst
                    permOrders{isub}(iperm,:) = [1 1+randperm(minSessions-1)];%ONLY PERMUTING FIRST 30 SESSIONS
                else
                    permOrders{isub}(iperm,:) = randperm(minSessions);%ONLY PERMUTING FIRST 30 SESSIONS
                end
            end
        end
        temp = squeeze(avgImgCorrMat(iroi,isub,1:minSessions,1:minSessions));
        tempOri = squeeze(avgImgCorrMatOri(iroi,isub,1:minSessions,1:minSessions));
        for idist=1:minSessions-1
            betweenSessImg(iroi,isub,idist) = nanmean(temp(distMatrix==idist));
            betweenSessImgOri(iroi,isub,idist) = nanmean(tempOri(distMatrix==idist));
        end
        if toNormalize
            initIndex = 1;
            betweenSessImg(iroi,isub,:) = (betweenSessImg(iroi,isub,:) - betweenSessImg(iroi,isub,initIndex))./abs(betweenSessImg(iroi,isub,:) + betweenSessImg(iroi,isub,initIndex));
            betweenSessImgOri(iroi,isub,:) = (betweenSessImgOri(iroi,isub,:) - betweenSessImgOri(iroi,isub,initIndex))./abs(betweenSessImgOri(iroi,isub,:) + betweenSessImgOri(iroi,isub,initIndex));
        end
        for iperm=1:nperms
            permOrder = permOrders{isub}(iperm,:);%different for each subject, but same for all ROIs
            tempPerm = squeeze(betweenSessCorr(iroi,isub,permOrder,permOrder));
            tempOriPerm = squeeze(betweenSessCorrOri(iroi,isub,permOrder,permOrder));
            for idist=1:minSessions-1
                betweenSessImgPerm(iroi,isub,iperm,idist) = nanmean(tempPerm(distMatrix==idist));
                betweenSessImgOriPerm(iroi,isub,iperm,idist) = nanmean(tempOriPerm(distMatrix==idist));
            end
            imgCorrMatPerm(iroi,isub,iperm,:,:) = tempPerm;
            imgCorrMatOriPerm(iroi,isub,iperm,:,:) = tempOriPerm;
        end
    end
    
    
    %correlations between RDMs
    for isub=1:length(subjects)
        temp = squeeze(betweenSessCorr(iroi,isub,1:minSessions,1:minSessions)); %includes NaNs for sessinos that don't exist
        tempOri = squeeze(betweenSessCorrOri(iroi,isub,1:minSessions,1:minSessions));%includes NaNs for sessinos that don't exist
        for idist=1:minSessions-1
            betweenSessDist(iroi,isub,idist) = nanmean(temp(distMatrix==idist));
            betweenSessDistOri(iroi,isub,idist) = nanmean(tempOri(distMatrix==idist));
        end
        if toNormalize
            initIndex = 1;
            betweenSessDist(iroi,isub,:) = (betweenSessDist(iroi,isub,:) - betweenSessDist(iroi,isub,initIndex))./abs(betweenSessDist(iroi,isub,:) + betweenSessDist(iroi,isub,initIndex));
            betweenSessDistOri(iroi,isub,:) = (betweenSessDistOri(iroi,isub,:) - betweenSessDistOri(iroi,isub,initIndex))./abs(betweenSessDistOri(iroi,isub,:) + betweenSessDistOri(iroi,isub,initIndex));
        end
        for iperm=1:nperms
            permOrder = permOrders{isub}(iperm,:);%different for each subject, but same for all ROIs
            tempPerm = squeeze(betweenSessCorr(iroi,isub,permOrder,permOrder));
            tempOriPerm = squeeze(betweenSessCorrOri(iroi,isub,permOrder,permOrder));
            for idist=1:nsessions-1
                betweenSessDistPerm(iroi,isub,iperm,idist) = nanmean(tempPerm(distMatrix==idist));
                betweenSessDistOriPerm(iroi,isub,iperm,idist) = nanmean(tempOriPerm(distMatrix==idist));
            end
            betweenSessCorrPerm(iroi,isub,iperm,:,:) = tempPerm;
            betweenSessCorrOriPerm(iroi,isub,iperm,:,:) = tempOriPerm;
        end
    end
end

save(fullfile(saveFolder, ['permPop' fixedFirstStr num2str(nperms) zscoreStr r2threshStr '.mat']),...
    'rois','toNormalize','toZscore','r2thresh','nrois', ...
    'permOrders', 'subSessions', 'subjects', 'minSessions','distMatrix',...
    'betweenSessCorr','betweenSessCorrOri',...
    'avgImgCorrMat','avgImgCorrMatOri',...
    'betweenSessImg','betweenSessImgOri','betweenSessImgPerm','betweenSessImgOriPerm',...
    'betweenSessDist','betweenSessDistOri','betweenSessDistPerm','betweenSessDistOriPerm',...
    'numGoodVox',...
    'imgCorrMatPerm','imgCorrMatOriPerm','betweenSessCorrPerm','betweenSessCorrOriPerm');
toc
