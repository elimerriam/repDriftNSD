% savePerms_expand.m
%
% associated with the following publication: Roth, ZN, and Merriam, EP (2023).
% Representations in human primary visual cortex drift over time
% DOI:
%
%   usage: savePerms_expand(1,0,1000,0)
%   by: zvi roth
%   date: 3/10/2022
%   purpose:  permute the session order and compute measures of representational drift
%   uses files created by: regressSessionCombineRoi_expand.m
%   creates files used by: 


function res = savePerms_expand(rois,toZscore,nperms,r2thresh,fixedFirst)

if ieNotDefined('rois'), rois = [1]; end
if ieNotDefined('toZscore'), toZscore = 0; end
if ieNotDefined('nperms'), nperms = 1000; end
if ieNotDefined('r2thresh'), r2thresh = 0; end
if ieNotDefined('fixedFirst'), fixedFirst = 0; end

fixedFirstStr='';
if fixedFirst
    fixedFirstStr = '_fixedFirst_';
end

tic
rng(1)

saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand','/');
if ~isfolder(saveFolder)
    saveFolder = ['/misc/data18/rothzn/nsd/repDrift_expand/'];
end

toNormalize = 0;% (r2distant-r2initial)/abs(r2distant+r2initial)

useMedian=1;

subjects = [1:8];

r2threshStr = '';
if r2thresh>0
   r2threshStr = ['r2thresh' num2str(r2thresh,'%4.0f')];
end

nsubjects=length(subjects);
nrois=4;
nsessions=30;
% maxSessions = nsessions;
minSessions=nsessions;

meanBetas = NaN(nrois,nsubjects,nsessions);
meanStdBetas = NaN(nrois,nsubjects,nsessions);
meanConstantCoef = NaN(nrois,nsubjects,nsessions);
meanConstantOriCoef = NaN(nrois,nsubjects,nsessions);

distMatrix = toeplitz(0:nsessions-1);%nsplits = nsessions+1. from 0 to 39.
lagMatrix = toeplitz(0:-1:1-nsessions,0:nsessions-1);%from -39 to +39
subSessions = zeros(nsubjects,nsessions);

rssDistSess = cell(nrois,nsessions);
rssOriDistSess = cell(nrois,nsessions);
r2DistSess = cell(nrois,nsessions);
r2OriDistSess = cell(nrois,nsessions);
pearsonDistSess = cell(nrois,nsessions);
pearsonOriDistSess = cell(nrois,nsessions);

rssSplit = cell(nrois,1);
rssOriSplit = cell(nrois,1);
r2split = cell(nrois,1);
r2oriSplit = cell(nrois,1);
pearsonRori = cell(nrois,1);
pearsonR = cell(nrois,1);

rssDist = cell(nrois,1);
rssOriDist = cell(nrois,1);
r2Dist = cell(nrois,1);
r2OriDist = cell(nrois,1);
pearsonDist = cell(nrois,1);
pearsonOriDist = cell(nrois,1);

rssDistPerm = cell(nrois,1);
rssOriDistPerm = cell(nrois,1);
r2DistPerm = cell(nrois,1);
r2OriDistPerm = cell(nrois,1);
pearsonDistPerm = cell(nrois,1);
pearsonOriDistPerm = cell(nrois,1);

for iroi=1:nrois
    rssSplit{iroi} = NaN(length(subjects),nsessions,nsessions);
    rssOriSplit{iroi} = NaN(length(subjects),nsessions,nsessions);
    r2split{iroi} = NaN(length(subjects),nsessions,nsessions);
    r2oriSplit{iroi} = NaN(length(subjects),nsessions,nsessions);
    pearsonRori{iroi} = NaN(length(subjects),nsessions,nsessions);
    pearsonR{iroi} = NaN(length(subjects),nsessions,nsessions);
    
    rssDist{iroi} = NaN(length(subjects),nsessions-1);
    rssOriDist{iroi} = NaN(length(subjects),nsessions-1);
    r2Dist{iroi} = NaN(length(subjects),nsessions-1);
    r2OriDist{iroi} = NaN(length(subjects),nsessions-1);
    pearsonDist{iroi} = NaN(length(subjects),nsessions-1);
    pearsonOriDist{iroi} = NaN(length(subjects),nsessions-1);
    
     rssDistPerm{iroi} = NaN(length(subjects),nperms,nsessions-1);
     rssOriDistPerm{iroi} = NaN(length(subjects),nperms,nsessions-1);
    r2DistPerm{iroi} = NaN(length(subjects),nperms,nsessions-1);
    r2OriDistPerm{iroi} = NaN(length(subjects),nperms,nsessions-1);
    pearsonDistPerm{iroi} = NaN(length(subjects),nperms,nsessions-1);
    pearsonOriDistPerm{iroi} = NaN(length(subjects),nperms,nsessions-1);
    
    rssPerm{iroi} = NaN(length(subjects),nperms,nsessions,nsessions);
    rssOriPerm{iroi} = NaN(length(subjects),nperms,nsessions,nsessions);
    r2perm{iroi} = NaN(length(subjects),nperms,nsessions,nsessions);
    r2oriPerm{iroi} = NaN(length(subjects),nperms,nsessions,nsessions);
    pearsonPerm{iroi} = NaN(length(subjects),nperms,nsessions,nsessions);
    pearsonOriPerm{iroi} = NaN(length(subjects),nperms,nsessions,nsessions);
    
end


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
autocorrLength = 21;
meanAutocorrBetas = zeros(nrois,nsubjects,autocorrLength);
autocorrMeanBetas = zeros(nrois,nsubjects,autocorrLength);
meanAutocorrStdBetas = zeros(nrois,nsubjects,autocorrLength);
autocorrMeanStdBetas = zeros(nrois,nsubjects,autocorrLength);

meanAutocorrConstant = zeros(nrois,nsubjects,autocorrLength);
meanAutocorrConstantOri = zeros(nrois,nsubjects,autocorrLength);
autocorrMeanConstant = zeros(nrois,nsubjects,autocorrLength);
autocorrMeanConstantOri = zeros(nrois,nsubjects,autocorrLength);

autocorrMeanCoef = zeros(nrois,nsubjects,autocorrLength);
autocorrMeanCoefOri = zeros(nrois,nsubjects,autocorrLength);
meanAutocorrCoef = zeros(nrois,nsubjects,autocorrLength);
meanAutocorrCoefOri = zeros(nrois,nsubjects,autocorrLength);

meanAutocorrBetasPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
meanAutocorrStdBetasPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
autocorrMeanBetasPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
autocorrMeanStdBetasPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
meanAutocorrConstantPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
meanAutocorrConstantOriPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
autocorrMeanConstantPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
autocorrMeanConstantOriPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
autocorrMeanCoefPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
autocorrMeanCoefOriPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
meanAutocorrCoefPerm = zeros(nrois,nsubjects,nperms,autocorrLength);
meanAutocorrCoefOriPerm = zeros(nrois,nsubjects,nperms,autocorrLength);

autocorrBetasPerm = cell(nrois,nsubjects);
autocorrConstantPerm = cell(nrois,nsubjects);
autocorrConstantOriPerm = cell(nrois,nsubjects);
autocorrCoefPerm = cell(nrois,nsubjects);
autocorrCoefOriPerm = cell(nrois,nsubjects);
load(fullfile(saveFolder,['betasVar.mat']),'roiBetasVar30');
roiBetasVar = roiBetasVar30;

for isub=1:nsubjects
    isub
    subjNum = subjects(isub);
    load(fullfile(saveFolder, ['regressSessCombineROI_sub' num2str(subjNum) zscoreStr '.mat']),'allRoiPrf',...
        'roiNsdCorr','roiNsdOriCorr','roiNsdOriR2','roiNsdR2',...
        'roiNsdOriRss','roiNsdRss',...
        'nsplits',...
        'sessBetas','sessStdBetas','voxConstCoef','voxConstOriCoef',...
        'voxOriCoef','voxCoef');

    
    nsessions = nsplits-1; %maxnumber of  sessions. nsplits includes the mean, and is <41 for some subjects
    subSessions(isub,1:nsessions) = ones;
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
    for iroi=rois
        %normalize RSS
        for isplit=1:nsplits
            roiNsdOriRss{iroi}(isplit,:,:) = squeeze(roiNsdOriRss{iroi}(isplit,:,:))./roiBetasVar{isub,iroi};
        end
        
        %for plotting visual field
        subRoiPrf{isub,iroi} = allRoiPrf{iroi};%for plotting visual field
        
        %if r2thresh>0 then use only voxels above threshold
        goodVox = allRoiPrf{iroi}.r2>r2thresh;
        r2GoodVox(iroi,isub) = r2thresh;
        totalVox(iroi,isub) = length(goodVox);
        numGoodVox(iroi,isub) = sum(goodVox);
        
        meanBetas(iroi,isub,1:nsessions) = mean(sessBetas{iroi}(1:nsessions,goodVox),2);
        meanStdBetas(iroi,isub,1:nsessions) = mean(sessStdBetas{iroi}(1:nsessions,goodVox),2);
        meanConstantCoef(iroi,isub,1:nsessions) = mean(voxConstCoef{iroi}(1:nsessions,goodVox),2);
        meanConstantOriCoef(iroi,isub,1:nsessions) = mean(voxConstOriCoef{iroi}(1:nsessions,goodVox),2);
        
        %std over coefficients excluding constant
        stdCoef{iroi,isub} = std(voxCoef{iroi}(:,:,1:end-1),0,3);%session,vox
        stdOriCoef{iroi,isub} = std(voxOriCoef{iroi}(:,:,1:end-1),0,3);
        
        if useMedian %median of voxels for each subject and ROI
            rssSplit{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(roiNsdRss{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            rssOriSplit{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(roiNsdOriRss{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            r2split{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(roiNsdR2{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            r2oriSplit{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(roiNsdOriR2{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            pearsonR{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(roiNsdCorr{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            pearsonRori{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(roiNsdOriCorr{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
        else %use mean
            rssSplit{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmean(roiNsdRss{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            rssOriSplit{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmean(roiNsdOriRss{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            r2split{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmean(roiNsdR2{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            r2oriSplit{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmean(roiNsdOriR2{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            pearsonR{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmean(roiNsdCorr{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
            pearsonRori{iroi}(isub,1:nsessions,1:nsessions) = squeeze(nanmean(roiNsdOriCorr{iroi}(1:nsessions,goodVox,:),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
        end

        for idist=1:nsessions
           temp = squeeze(rssSplit{iroi}(isub,1:minSessions,1:minSessions));
            rssDistSess{iroi,idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
            rssDist{iroi}(isub,idist) = nanmean(temp(distMatrix==idist));

            temp = squeeze(rssOriSplit{iroi}(isub,1:minSessions,1:minSessions));
            rssOriDistSess{iroi,idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
            rssOriDist{iroi}(isub,idist) = nanmean(temp(distMatrix==idist));

            temp = squeeze(r2split{iroi}(isub,1:minSessions,1:minSessions));%includes NaNs for sessinos that don't exist
            r2DistSess{iroi,idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
            r2Dist{iroi}(isub,idist) = nanmean(temp(distMatrix==idist));
            
            temp = squeeze(r2oriSplit{iroi}(isub,1:minSessions,1:minSessions));
            r2OriDistSess{iroi,idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
            r2OriDist{iroi}(isub,idist) = nanmean(temp(distMatrix==idist));
            
            temp = squeeze(pearsonR{iroi}(isub,1:minSessions,1:minSessions));
            pearsonDistSess{iroi,idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
            pearsonDist{iroi}(isub,idist) = nanmean(temp(distMatrix==idist));
            
            temp = squeeze(pearsonRori{iroi}(isub,1:minSessions,1:minSessions));
            pearsonOriDistSess{iroi,idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
            pearsonOriDist{iroi}(isub,idist) = nanmean(temp(distMatrix==idist));
            
        end
        
        if toNormalize
%             initIndex = 2;%first index is within session
            initIndex = 1;%first index is lag of one
            rssOriDist{iroi}(isub,:) = (rssOriDist{iroi}(isub,:) - rssOriDist{iroi}(isub,initIndex))./abs(rssOriDist{iroi}(isub,:) + rssOriDist{iroi}(isub,initIndex));
            rssDist{iroi}(isub,:) = (rssDist{iroi}(isub,:) - rssDist{iroi}(isub,initIndex))./abs(rssDist{iroi}(isub,:) + rssDist{iroi}(isub,initIndex));
            r2Dist{iroi}(isub,:) = (r2Dist{iroi}(isub,:) - r2Dist{iroi}(isub,initIndex))./abs(r2Dist{iroi}(isub,:) + r2Dist{iroi}(isub,initIndex));
            r2OriDist{iroi}(isub,:) = (r2OriDist{iroi}(isub,:) - r2OriDist{iroi}(isub,initIndex))./abs(r2OriDist{iroi}(isub,:) + r2OriDist{iroi}(isub,initIndex));
            pearsonDist{iroi}(isub,:) = (pearsonDist{iroi}(isub,:) - pearsonDist{iroi}(isub,initIndex))./abs(pearsonDist{iroi}(isub,:) + pearsonDist{iroi}(isub,initIndex));
            pearsonOriDist{iroi}(isub,:) = (pearsonOriDist{iroi}(isub,:) - pearsonOriDist{iroi}(isub,initIndex))./abs(pearsonOriDist{iroi}(isub,:) + pearsonOriDist{iroi}(isub,initIndex));
        end
        
        for ivox=1:size(sessBetas{iroi},2)
            autocorrBetas{iroi,isub}(ivox,:) = autocorr(sessBetas{iroi}(:,ivox));
            autocorrStdBetas{iroi,isub}(ivox,:) = autocorr(sessStdBetas{iroi}(:,ivox));
            autocorrConstant{iroi,isub}(ivox,:) = autocorr(voxConstCoef{iroi}(:,ivox));
            autocorrConstantOri{iroi,isub}(ivox,:) = autocorr(voxConstOriCoef{iroi}(:,ivox));
            autocorrCoef{iroi,isub}(ivox,:) = autocorr(stdCoef{iroi,isub}(:,ivox));
            autocorrCoefOri{iroi,isub}(ivox,:) = autocorr(stdOriCoef{iroi,isub}(:,ivox));
        end
        meanAutocorrBetas(iroi,isub,:) = squeeze(mean(autocorrBetas{iroi,isub}(goodVox,:),1));
        autocorrMeanBetas(iroi,isub,:) = autocorr(squeeze(meanBetas(iroi,isub,1:nsessions)));
        meanAutocorrStdBetas(iroi,isub,:) = squeeze(mean(autocorrStdBetas{iroi,isub}(goodVox,:),1));
        autocorrMeanStdBetas(iroi,isub,:) = autocorr(squeeze(meanStdBetas(iroi,isub,1:nsessions)));
        
        meanAutocorrConstant(iroi,isub,:) = squeeze(mean(autocorrConstant{iroi,isub}(goodVox,:),1));
        meanAutocorrConstantOri(iroi,isub,:) = squeeze(mean(autocorrConstantOri{iroi,isub}(goodVox,:),1));
        autocorrMeanConstant(iroi,isub,:) = autocorr(squeeze(meanConstantCoef(iroi,isub,1:nsessions)));
        autocorrMeanConstantOri(iroi,isub,:) = autocorr(squeeze(meanConstantOriCoef(iroi,isub,1:nsessions)));
        
        autocorrMeanCoef(iroi,isub,:) = autocorr(squeeze(mean(stdCoef{iroi,isub}(:,goodVox),2)));%mean over voxels
        autocorrMeanCoefOri(iroi,isub,:) = autocorr(squeeze(mean(stdOriCoef{iroi,isub}(:,goodVox),2)));%mean over voxels
        meanAutocorrCoef(iroi,isub,:) = squeeze(nanmean(autocorrCoef{iroi,isub}(goodVox,:),1));
        meanAutocorrCoefOri(iroi,isub,:) = squeeze(nanmean(autocorrCoefOri{iroi,isub}(goodVox,:),1));

        autocorrBetasPerm{iroi,isub} = zeros(nperms,size(sessBetas{iroi},2),autocorrLength);
        autocorrConstantPerm{iroi,isub} = zeros(nperms,size(sessBetas{iroi},2),autocorrLength);
        autocorrConstantOriPerm{iroi,isub} = zeros(nperms,size(sessBetas{iroi},2),autocorrLength); 
        autocorrCoefPerm{iroi,isub} = zeros(nperms,size(sessBetas{iroi},2),autocorrLength); 
        autocorrCoefOriPerm{iroi,isub} = zeros(nperms,size(sessBetas{iroi},2),autocorrLength);
        
        
        for iperm=1:nperms
            permOrder = permOrders{isub}(iperm,:);%different for each subject, but same for all ROIs

            for idist=1:nsessions-1
                temp = squeeze(rssSplit{iroi}(isub,permOrder,permOrder));
                rssDistSessPerm{iroi,idist}(isub,iperm,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
                rssDistPerm{iroi}(isub,iperm,idist) = nanmean(temp(distMatrix==idist));
                rssPerm{iroi}(isub,iperm,:,:) = temp;

                temp = squeeze(rssOriSplit{iroi}(isub,permOrder,permOrder));
                rssOriDistSessPerm{iroi,idist}(isub,iperm,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
                rssOriDistPerm{iroi}(isub,iperm,idist) = nanmean(temp(distMatrix==idist));
                rssOriPerm{iroi}(isub,iperm,:,:) = temp;

                temp = squeeze(r2split{iroi}(isub,permOrder,permOrder));
                r2DistSessPerm{iroi,idist}(isub,iperm,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
                r2DistPerm{iroi}(isub,iperm,idist) = nanmean(temp(distMatrix==idist));
                r2perm{iroi}(isub,iperm,:,:) = temp;
                
                temp = squeeze(r2oriSplit{iroi}(isub,permOrder,permOrder));
                r2OriDistSessPerm{iroi,idist}(isub,iperm,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
                r2OriDistPerm{iroi}(isub,iperm,idist) = nanmean(temp(distMatrix==idist));
                r2oriPerm{iroi}(isub,iperm,:,:) = temp;
                
                temp = squeeze(pearsonR{iroi}(isub,permOrder,permOrder));
                pearsonDistSessPerm{iroi,idist}(isub,iperm,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
                pearsonDistPerm{iroi}(isub,iperm,idist) = nanmean(temp(distMatrix==idist));
                pearsonPerm{iroi}(isub,iperm,:,:) = temp;
                
                temp = squeeze(pearsonRori{iroi}(isub,permOrder,permOrder));
                pearsonOriDistSessPerm{iroi,idist}(isub,iperm,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
                pearsonOriDistPerm{iroi}(isub,iperm,idist) = nanmean(temp(distMatrix==idist));
                pearsonOriPerm{iroi}(isub,iperm,:,:) = temp;
            end

            if toNormalize
                initIndex = 1;
                rssDistPerm{iroi}(isub,iperm,:) = (rssDistPerm{iroi}(isub,iperm,:) - rssDistPerm{iroi}(isub,iperm,initIndex))./abs(rssDistPerm{iroi}(isub,iperm,:) + rssDistPerm{iroi}(isub,iperm,initIndex));
                rssOriDistPerm{iroi}(isub,iperm,:) = (rssOriDistPerm{iroi}(isub,iperm,:) - rssOriDistPerm{iroi}(isub,iperm,initIndex))./abs(rssOriDistPerm{iroi}(isub,iperm,:) + rssOriDistPerm{iroi}(isub,iperm,initIndex));
                r2DistPerm{iroi}(isub,iperm,:) = (r2DistPerm{iroi}(isub,iperm,:) - r2DistPerm{iroi}(isub,iperm,initIndex))./abs(r2DistPerm{iroi}(isub,iperm,:) + r2DistPerm{iroi}(isub,iperm,initIndex));
                r2OriDistPerm{iroi}(isub,iperm,:) = (r2OriDistPerm{iroi}(isub,iperm,:) - r2OriDistPerm{iroi}(isub,iperm,initIndex))./abs(r2OriDistPerm{iroi}(isub,iperm,:) + r2OriDistPerm{iroi}(isub,iperm,initIndex));
                pearsonDistPerm{iroi}(isub,iperm,:) = (pearsonDistPerm{iroi}(isub,iperm,:) - pearsonDistPerm{iroi}(isub,iperm,initIndex))./abs(pearsonDistPerm{iroi}(isub,iperm,:) + pearsonDistPerm{iroi}(isub,iperm,initIndex));
                pearsonOriDistPerm{iroi}(isub,iperm,:) = (pearsonOriDistPerm{iroi}(isub,iperm,:) - pearsonOriDistPerm{iroi}(isub,iperm,initIndex))./abs(pearsonOriDistPerm{iroi}(isub,iperm,:) + pearsonOriDistPerm{iroi}(isub,iperm,initIndex));
            end
            for ivox=1:size(sessBetas{iroi},2)
                autocorrBetasPerm{iroi,isub}(iperm,ivox,:) = autocorr(sessBetas{iroi}(permOrder,ivox));
                autocorrStdBetasPerm{iroi,isub}(iperm,ivox,:) = autocorr(sessStdBetas{iroi}(permOrder,ivox));
                autocorrConstantPerm{iroi,isub}(iperm,ivox,:) = autocorr(voxConstCoef{iroi}(permOrder,ivox));
                autocorrConstantOriPerm{iroi,isub}(iperm,ivox,:) = autocorr(voxConstOriCoef{iroi}(permOrder,ivox));
                autocorrCoefPerm{iroi,isub}(iperm,ivox,:) = autocorr(stdCoef{iroi,isub}(permOrder,ivox));
                autocorrCoefOriPerm{iroi,isub}(iperm,ivox,:) = autocorr(stdOriCoef{iroi,isub}(permOrder,ivox));
            end
            meanAutocorrBetasPerm(iroi,isub,iperm,:) = squeeze(mean(autocorrBetasPerm{iroi,isub}(iperm,goodVox,:),2));
            meanAutocorrStdBetasPerm(iroi,isub,iperm,:) = squeeze(mean(autocorrStdBetasPerm{iroi,isub}(iperm,goodVox,:),2));
            autocorrMeanBetasPerm(iroi,isub,iperm,:) = autocorr(squeeze(meanBetas(iroi,isub,permOrder)));
            autocorrMeanStdBetasPerm(iroi,isub,iperm,:) = autocorr(squeeze(meanStdBetas(iroi,isub,permOrder)));

            meanAutocorrConstantPerm(iroi,isub,iperm,:) = squeeze(mean(autocorrConstantPerm{iroi,isub}(iperm,goodVox,:),2));
            meanAutocorrConstantOriPerm(iroi,isub,iperm,:) = squeeze(mean(autocorrConstantOriPerm{iroi,isub}(iperm,goodVox,:),2));
            autocorrMeanConstantPerm(iroi,isub,iperm,:) = autocorr(squeeze(meanConstantCoef(iroi,isub,permOrder)));
            autocorrMeanConstantOriPerm(iroi,isub,iperm,:) = autocorr(squeeze(meanConstantOriCoef(iroi,isub,permOrder)));
            
            autocorrMeanCoefPerm(iroi,isub,iperm,:) = autocorr(squeeze(mean(stdCoef{iroi,isub}(permOrder,goodVox),2)));%mean over voxels
            autocorrMeanCoefOriPerm(iroi,isub,iperm,:) = autocorr(squeeze(mean(stdOriCoef{iroi,isub}(permOrder,goodVox),2)));%mean over voxels
            meanAutocorrCoefPerm(iroi,isub,iperm,:) = squeeze(nanmean(autocorrCoefPerm{iroi,isub}(iperm,goodVox,:),2));
            meanAutocorrCoefOriPerm(iroi,isub,iperm,:) = squeeze(nanmean(autocorrCoefOriPerm{iroi,isub}(iperm,goodVox,:),2));
            
        end
        
    end
    
end

normalizeStr = '';
if toNormalize
    normalizeStr = '_normalized';
end

save(fullfile(saveFolder, ['perm' fixedFirstStr num2str(nperms) zscoreStr normalizeStr r2threshStr  '.mat']),...
    'toNormalize','toZscore', 'useMedian','r2thresh','nrois','rois', ...
    'permOrders', 'subSessions', 'subjects', 'minSessions',...
    'r2split','r2oriSplit','pearsonRori','pearsonR',...
    'r2Dist','r2OriDist','pearsonDist','pearsonOriDist',...
    'r2DistPerm','r2OriDistPerm','pearsonDistPerm','pearsonOriDistPerm',....
    'r2DistSess','r2OriDistSess','pearsonDistSess','pearsonOriDistSess',...
    'rssDistSess','rssOriDistSess',...
    'rssDist','rssOriDist',...
    'rssSplit','rssOriSplit',...
    'rssDistPerm','rssOriDistPerm',...
    'rssDistSessPerm','rssOriDistSessPerm',...
    'r2DistSessPerm','r2OriDistSessPerm','pearsonDistSessPerm','pearsonOriDistSessPerm',...
    'meanAutocorrBetas', 'autocorrMeanBetas','meanAutocorrStdBetas','autocorrMeanStdBetas',...
    'meanAutocorrConstant','meanAutocorrConstantOri','autocorrMeanConstant','autocorrMeanConstantOri',...
    'autocorrMeanCoef','autocorrMeanCoefOri','meanAutocorrCoef','meanAutocorrCoefOri',...
    'meanAutocorrBetasPerm','autocorrMeanBetasPerm','meanBetas',...
    'meanAutocorrStdBetasPerm','autocorrMeanStdBetasPerm',...
    'meanAutocorrConstantPerm','meanAutocorrConstantOriPerm',...
    'autocorrMeanConstantPerm','autocorrMeanConstantOriPerm',...
    'autocorrMeanCoefPerm','autocorrMeanCoefOriPerm',...
    'meanAutocorrCoefPerm','meanAutocorrCoefOriPerm',...
    'subRoiPrf','numGoodVox','totalVox',...
    'r2perm','r2oriPerm','pearsonPerm','pearsonOriPerm',...
    'rssPerm','rssOriPerm');

toc

res = fullfile(saveFolder, ['perm' num2str(nperms) zscoreStr normalizeStr r2threshStr  '.mat']);