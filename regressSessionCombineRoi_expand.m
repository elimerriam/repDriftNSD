% regressSessionCombineRoi_expand.m
%
% associated with the following publication: Roth, ZN, and Merriam, EP (2023).
% Representations in human primary visual cortex drift over time
% DOI:
%
%   usage: for isubject=1:8; regressSessionCombineRoi_expand(isubject,1,0); end;
%   by: zvi roth
%   date: 3/10/2022
%   purpose:  combine data from ventral and dorsal ROIs
%   uses files created by: regressPrfSplit_expand.m
%   creates files used by: savePerms_expand.m; simPopResponse_expand

function res = regressSessionCombineRoi_expand(isub,numregions,toZscore)
if ieNotDefined('numregions'), numregions = 1; end
if ieNotDefined('toZscore'), toZscore = 0; end

saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand');
if ~isfolder(saveFolder)
    saveFolder = '/misc/data18/rothzn/nsd/repDrift_expand/';
end

if ieNotDefined('toZscore'), toZscore = 0; end
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
for iregion=1:numregions
    
    visualRegion = iregion;%V1,V2,V3,V4
    
    load(fullfile(saveFolder,['regressPrfSplit_session_v' num2str(visualRegion) '_sub' num2str(isub) zscoreStr  '.mat']), ...
        'nsd', ... %'synth',...
        'rois','roiPrf','nsplits');
    if length(rois)>1
        oldNsd = nsd;
        clear nsd
        
        nsd.sessBetas{1} = [];
        nsd.sessStdBetas{1} = [];
        
        nsd.voxOriCoefCorrWithConst{1} = [];
        nsd.voxCoefCorrWithConst{1} = [];
        nsd.voxOriCoefCorr{1} = [];
        nsd.voxCoefCorr{1} = [];
    
        nsd.r2{1} = [];
        nsd.r2ori{1} = [];
        nsd.r2split{1} = [];
        nsd.r2oriSplit{1} = [];
        nsd.rssSplit{1} = [];
        nsd.rssOriSplit{1} = [];
        nsd.voxCoef{1} = [];
        nsd.voxOriCoef{1} = [];
        nsd.voxPredOriCoef{1} = [];
        nsd.voxOriPredOriCoef{1} = [];
        nsd.voxResidOriCoef{1} = [];
        nsd.voxOriResidOriCoef{1} = [];
        
        nsd.voxPredOriR2{1} = [];
        nsd.voxOriPredOriR2{1} = [];
        nsd.voxResidOriR2{1} = [];
        nsd.voxOriResidOriR2{1} = [];
        
        nsd.pearsonRori{1} = [];
        nsd.pearsonR{1} = [];
        
        nsd.roiInd{1} = [];

        
        for iroi=1:length(rois)%rois=2
            
            nsd.sessBetas{1} = cat(2,nsd.sessBetas{1}, oldNsd.sessBetas{iroi});
            nsd.sessStdBetas{1} = cat(2,nsd.sessStdBetas{1}, oldNsd.sessStdBetas{iroi});
            
            nsd.roiInd{1} = cat(1,nsd.roiInd{1}, oldNsd.roiInd{iroi});
            nsd.pearsonRori{1} = cat(2,nsd.pearsonRori{1},oldNsd.pearsonRori{iroi});
            nsd.pearsonR{1} = cat(2,nsd.pearsonR{1},oldNsd.pearsonR{iroi});
            nsd.r2split{1} = cat(2,nsd.r2split{1},oldNsd.r2split{iroi});
            nsd.r2oriSplit{1} = cat(2,nsd.r2oriSplit{1},oldNsd.r2oriSplit{iroi});
            nsd.rssSplit{1} = cat(2,nsd.rssSplit{1},oldNsd.rssSplit{iroi});
            nsd.rssOriSplit{1} = cat(2,nsd.rssOriSplit{1},oldNsd.rssOriSplit{iroi});
            nsd.voxCoef{1} = cat(2,nsd.voxCoef{1},oldNsd.voxCoef{iroi});
            nsd.voxOriCoef{1} = cat(2,nsd.voxOriCoef{1},oldNsd.voxOriCoef{iroi});
            nsd.voxPredOriCoef{1} = cat(2,nsd.voxPredOriCoef{1},oldNsd.voxPredOriCoef{iroi});
            nsd.voxOriPredOriCoef{1} = cat(2,nsd.voxOriPredOriCoef{1},oldNsd.voxOriPredOriCoef{iroi});
            nsd.voxResidOriCoef{1} = cat(2,nsd.voxResidOriCoef{1},oldNsd.voxResidOriCoef{iroi});
            nsd.voxOriResidOriCoef{1} = cat(2,nsd.voxOriResidOriCoef{1},oldNsd.voxOriResidOriCoef{iroi});
            
            nsd.r2{1} = cat(2,nsd.r2{1},oldNsd.r2{iroi});
            nsd.r2ori{1} = cat(2,nsd.r2ori{1},oldNsd.r2ori{iroi});
            
            nsd.voxOriCoefCorrWithConst{1} = cat(1,nsd.voxOriCoefCorrWithConst{1},oldNsd.voxOriCoefCorrWithConst{iroi});
            nsd.voxCoefCorrWithConst{1} = cat(1,nsd.voxCoefCorrWithConst{1},oldNsd.voxCoefCorrWithConst{iroi});
            nsd.voxOriCoefCorr{1} = cat(1,nsd.voxOriCoefCorr{1},oldNsd.voxOriCoefCorr{iroi});
            nsd.voxCoefCorr{1} = cat(1,nsd.voxCoefCorr{1},oldNsd.voxCoefCorr{iroi});
        
            nsd.voxPredOriR2{1} = cat(2,nsd.voxPredOriR2{1},oldNsd.voxPredOriR2{iroi});
            nsd.voxOriPredOriR2{1} = cat(2,nsd.voxOriPredOriR2{1},oldNsd.voxOriPredOriR2{iroi});
            nsd.voxResidOriR2{1} = cat(2,nsd.voxResidOriR2{1},oldNsd.voxResidOriR2{iroi});
            nsd.voxOriResidOriR2{1} = cat(2,nsd.voxOriResidOriR2{1},oldNsd.voxOriResidOriR2{iroi});
 
        end
        oldPrf = roiPrf; clear roiPrf;
        roiPrf{1}.ecc=[];
        roiPrf{1}.ang=[];
        roiPrf{1}.sz=[];
        roiPrf{1}.r2=[];
        roiPrf{1}.x=[];
        roiPrf{1}.y=[];
        for iroi=1:length(rois)
            roiPrf{1}.ecc = cat(1,roiPrf{1}.ecc,oldPrf{iroi}.ecc);
            roiPrf{1}.ang = cat(1,roiPrf{1}.ang,oldPrf{iroi}.ang);
            roiPrf{1}.sz = cat(1,roiPrf{1}.sz,oldPrf{iroi}.sz);
            roiPrf{1}.r2 = cat(1,roiPrf{1}.r2,oldPrf{iroi}.r2);
            roiPrf{1}.x = cat(1,roiPrf{1}.x,oldPrf{iroi}.x);
            roiPrf{1}.y = cat(1,roiPrf{1}.y,oldPrf{iroi}.y);
        end
        rois = 1;
    end
    iroi=1;
    
    %% AVERAGE SPLITS
    
    nsd.voxOriCoefCorrWithConst{1}(nsplits+1,:,:) = mean(nsd.voxOriCoefCorrWithConst{1},1);
    nsd.voxCoefCorrWithConst{1}(nsplits+1,:,:) = mean(nsd.voxCoefCorrWithConst{1},1);
    nsd.voxOriCoefCorr{1}(nsplits+1,:,:) = mean(nsd.voxOriCoefCorr{1},1);
    nsd.voxCoefCorr{1}(nsplits+1,:,:) = mean(nsd.voxCoefCorr{1},1);
            
    nsd.voxCoef{1}(nsplits+1,:,:) = mean(nsd.voxCoef{1},1);
    nsd.voxOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriCoef{1},1);
    nsd.voxPredOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxPredOriCoef{1},1);
    nsd.voxOriPredOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriPredOriCoef{1},1);
    nsd.voxResidOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxResidOriCoef{1},1);
    nsd.voxOriResidOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriResidOriCoef{1},1);
    
    nsd.pearsonRori{1}(nsplits+1,:,:) = mean(nsd.pearsonRori{1},1);
    nsd.pearsonR{1}(nsplits+1,:,:) = mean(nsd.pearsonR{1},1);
    nsd.r2{1}(nsplits+1,:) = mean(nsd.r2{1},1);
    nsd.r2ori{1}(nsplits+1,:) = mean(nsd.r2ori{1},1);
    nsd.r2split{1}(nsplits+1,:,:) = mean(nsd.r2split{1},1);
    nsd.r2oriSplit{1}(nsplits+1,:,:) = mean(nsd.r2oriSplit{1},1);
    nsd.rssSplit{1}(nsplits+1,:,:) = mean(nsd.rssSplit{1},1);
    nsd.rssOriSplit{1}(nsplits+1,:,:) = mean(nsd.rssOriSplit{1},1);

    nsplits = nsplits+1;
    

    %save model coefficients
    voxOriCoef{iregion} = nsd.voxOriCoef{iroi};
    voxCoef{iregion} = nsd.voxCoef{iroi};
    
    allRoiPrf{iregion} = roiPrf{iroi};
    roiInd{iregion} = nsd.roiInd{iroi};
    roiNsdCorr{iregion} = nsd.pearsonR{iroi};
    roiNsdOriCorr{iregion} = nsd.pearsonRori{iroi};
    roiNsdOriR2{iregion} = nsd.r2oriSplit{iroi};
    roiNsdR2{iregion} = nsd.r2split{iroi};
    roiNsdOriRss{iregion} = nsd.rssOriSplit{iroi};
    roiNsdRss{iregion} = nsd.rssSplit{iroi};
    roiNsdR2within{iregion} = nsd.r2{iroi};
    roiNsdOriR2within{iregion} = nsd.r2ori{iroi};
    
    %constant conefficent
    voxConstCoef{iregion} = squeeze(nsd.voxCoef{iroi}(:,:,end));
    voxConstOriCoef{iregion} = squeeze(nsd.voxOriCoef{iroi}(:,:,end));
    
    voxOriCoefCorrWithConst{iregion} = nsd.voxOriCoefCorrWithConst{iroi};
    voxCoefCorrWithConst{iregion} = nsd.voxCoefCorrWithConst{iroi};
    voxOriCoefCorr{iregion} = nsd.voxOriCoefCorr{iroi};
    voxCoefCorr{iregion} = nsd.voxCoefCorr{iroi};
    
    sessBetas{iregion} = nsd.sessBetas{iroi};
    sessStdBetas{iregion} = nsd.sessStdBetas{iroi};
    
    roiNsdOriPredOriR2{iregion} = nsd.voxOriPredOriR2{iroi};
    roiNsdOriResidOriR2{iregion} = nsd.voxOriResidOriR2{iroi};
    roiNsdPredOriR2{iregion} = nsd.voxPredOriR2{iroi};
    roiNsdResidOriR2{iregion} = nsd.voxResidOriR2{iroi};

end

    save(fullfile(saveFolder, ['regressSessCombineROI_sub' num2str(isub) zscoreStr '.mat']),'allRoiPrf',...
        'roiNsdCorr','roiNsdOriCorr','roiNsdOriR2','roiNsdR2',...
        'roiNsdOriRss','roiNsdRss','numregions',...
        'roiNsdResidOriR2','roiNsdOriResidOriR2','roiNsdPredOriR2','roiNsdOriPredOriR2',...
        'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis','nsplits',...
        'roiNsdOriR2within','roiNsdR2within',...
        'voxOriCoefCorrWithConst', 'voxCoefCorrWithConst','voxOriCoefCorr','voxCoefCorr',...
        'sessBetas','sessStdBetas','voxConstCoef','voxConstOriCoef',...
        'voxOriCoef','voxCoef');

end
