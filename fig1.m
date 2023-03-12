% fig1.m
%
% associated with the following publication: Roth, ZN, and Merriam, EP (2023).
% Representations in human primary visual cortex drift over time
% DOI:
%
%   usage: fig1()
%   by: zvi roth
%   date: 3/10/2022
%   purpose:  Plots for Fig 1 and Fig 2
%   uses files created by: savePerms_expand.m
%   creates pdf files 

close all
clear all
%uses data saved by savePerms_expand.m
tic

saveFigs=0;

addColorbars = 0;
fixedFirst=0;

%For Fig 2 choose toZscore=2, and toZscore=3;
toZscore=0;%0=none, 1=zscore, 2=zero mean, 3=normalized std, 4=zero ROI mean

r2thresh = 0;

singleSubject=1;

coreyVersion = 2;
addScatter = 0;
subjects = [1:8];

nperms=1000;
histBins = 20;
scatterSize=5;
scatterFill = 'filled';

fixedFirstStr='';
if fixedFirst
    fixedFirstStr = '_fixedFirst_';
end


toNormalize = 0;

figRoi=1;
linewidthWide=2;
linewidthNarrow = 1;
errorbarColor = [0.2 0.2 0.2];
surfaceAlpha = 0.1;
subColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840],...
    [0.4660 0.6740 0.1880]/2+[0.8500 0.3250 0.0980]/2};
graySubColor = 0.7*ones(1,3);
%histogram color
faceColor = graySubColor;
edgeColor = 'none';

rois = [1];
versionStr = '';
if version==2
    versionStr = '2';
end
saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand');
figsFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand_figs');
if ~isfolder(saveFolder)
    saveFolder = ['/misc/data18/rothzn/nsd/repDrift_expand/'];
    figsFolder = ['/misc/data18/rothzn/nsd/repDrift_expand_figs/'];
end

colormapName = 'parula';
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


r2threshStr = '';
if r2thresh>0
    r2threshStr = ['r2thresh' num2str(r2thresh,'%4.0f')];
end
load(fullfile(saveFolder, ['perm' fixedFirstStr num2str(nperms) zscoreStr  r2threshStr '.mat']),...
    'toNormalize','toZscore', 'useMedian','r2thresh','nrois','rois', ...
    'permOrders', 'subSessions',...
    'r2split','r2oriSplit','pearsonRori','pearsonR',...
    'r2Dist','r2OriDist','pearsonDist','pearsonOriDist',...
    'r2DistPerm','r2OriDistPerm','pearsonDistPerm','pearsonOriDistPerm',....
    'r2DistSess','r2OriDistSess','pearsonDistSess','pearsonOriDistSess',...
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
    'subRoiPrf','numGoodVox',...
    'r2perm','r2oriPerm','pearsonPerm','pearsonOriPerm');
nsubjects = length(subjects);
ifig=0;

%% use data for selected subjects
subSessions = subSessions(subjects,:);
minSessions = min(sum(subSessions,2));

%% distance matrix
sess = 1:minSessions;
sessDiff = sess - sess';
nsessions = 40;
distMatrix = toeplitz(0:nsessions-1);
sessDistVec = distMatrix(1:minSessions,1:minSessions);
%%
ifig=ifig+1;

iroi=figRoi;%nrois
f=figure(ifig); ifig=ifig+1;  isubplot=0;
rows=4;
cols=3;

for i=1:2
    switch i
        case 1
            similarityData = r2oriSplit{iroi};
            distSessData = r2OriDistSess;
            distData = r2OriDist{iroi};
            titleStr = 'R^2';
            ylabelStr = 'cvR^2';
            distSessDataPerm = r2OriDistSessPerm;
            distDataPerm = r2OriDistPerm{iroi};
            similarityPerm = r2oriPerm{iroi};
        case 2
            similarityData = pearsonRori{iroi};
            distSessData = pearsonOriDistSess;
            distData = pearsonOriDist{iroi};
            titleStr = 'Pearson''s r';
            ylabelStr = 'Pearson''s r';
            distSessDataPerm = pearsonOriDistSessPerm;
            distDataPerm = pearsonOriDistPerm{iroi};
            similarityPerm = pearsonOriPerm{iroi};
            
    end
    %keep data for chosen subjects only
    similarityData = similarityData(subjects,:,:);
    distData = distData(subjects,:);
    distDataPerm = distDataPerm(subjects,:,:);
    %keep minimum sessions that all subjects have
    similarityData = similarityData(:,1:minSessions,1:minSessions);
    distSessData = distSessData(:,1:minSessions);
    distData = distData(:,1:minSessions-1);
    distSessDataPerm = distSessDataPerm{1}(subjects,:,1:minSessions-1);
    distDataPerm = distDataPerm(:,:,1:minSessions-1);
    
    %matrix: trained on each session and tested on every session
    subplot(rows,cols,1+(i-1)*2*cols);
    temp = squeeze(mean(similarityData));
    
    img=imagesc(temp,'AlphaData',abs(sessDiff)>0); axis square%mean across subjects
    caxis([min(temp(abs(sessDiff)>0)), max(temp(abs(sessDiff)>0))]);

    if ~saveFigs
        title(titleStr);
    end
    colormap(colormapName);
    if addColorbars
        colorbar('southoutside')
        xticks([]);
        yticks([]);
    else
        
        ylabel('train session');
        xlabel('test session');
        xticks([1 minSessions]);
        yticks([1 minSessions]);
    end
    
    %similarity as function of distance, averaged over all possible
    %train-test pairs
    subplot(rows,cols,2+(i-1)*2*cols); axis square; hold all
    
    switch coreyVersion
        case 0
            [distCorr(i), pDistCorr(i)] = corr([1:minSessions-1]',mean(distData)');
            for iperm=1:nperms
                [distCorrPerm(i,iperm), pDistCorrPerm(i,iperm)] = corr([1:minSessions-1]',mean(squeeze(distDataPerm(:,iperm,:)))');
            end
            pPermDistCorr(i) = sum(distCorrPerm(i,:)<=distCorr(i))/nperms;
            %plot
            if singleSubject
                for isub=1:length(subjects)
                    plot(distData(isub,:),'color',subColor{subjects(isub)},'linewidth',linewidthNarrow); hold all
                end
            else
                dsErrorsurface(1:size(distData,2),mean(distData),std(distData)/sqrt(nsubjects),errorbarColor,surfaceAlpha);
            end
            plot(mean(distData),'linewidth', linewidthWide,'color','k'); hold all
            
        case 1
            %compute correlation using all elements of the goodness-of-fit matrix:
            dataDistVec = squeeze(mean(similarityData,1));%mean over subjects
            [distCorr(i), pDistCorr(i)] = corr(sessDistVec(sessDistVec>0),dataDistVec(sessDistVec>0));
            permDistMat = squeeze(mean(similarityPerm,1)); %mean over subjects, for each permutation
            for iperm=1:nperms
                permDistVec = permDistMat(iperm,1:minSessions,1:minSessions);
                [distCorrPerm(i,iperm), pDistCorrPerm(i,iperm)] = corr(sessDistVec(sessDistVec>0),permDistVec(sessDistVec>0));
            end
            pPermDistCorr(i) = sum(distCorrPerm(i,:)<=distCorr(i))/nperms;
            %plot
            scatter1=scatter(sessDistVec(sessDistVec(:)>0),dataDistVec(sessDistVec(:)>0),scatterSize,[0 0 0],scatterFill);
            scatter1.MarkerFaceAlpha = 2*surfaceAlpha;
            scatter1.MarkerEdgeAlpha = 2*surfaceAlpha;
            for idist=1:minSessions-1
                meanDistData(idist) = mean(dataDistVec(sessDistVec(:)==idist));
            end
            plot(meanDistData,'linewidth', linewidthWide,'color','k'); hold all
        case 2
            %compute correlation per subject and then average
            for isub=1:length(subjects)
                subData = similarityData(isub,1:minSessions,1:minSessions);
                [subDistCorr(isub,i) pSubDistCorr(isub,i)] = corr(sessDistVec(sessDistVec>0),subData(sessDistVec>0));
                permDistMat = squeeze(similarityPerm(isub,:,:,:));
                for iperm=1:nperms
                    permDistVec = permDistMat(iperm,1:minSessions,1:minSessions);
                    [subDistCorrPerm(i,isub,iperm), pSubDistCorrPerm(i,isub,iperm)] = corr(sessDistVec(sessDistVec>0),permDistVec(sessDistVec>0));
                end
                %plot
                if addScatter
                    scatter1 = scatter(sessDistVec(sessDistVec(:)>0),subData(sessDistVec(:)>0),scatterSize,repmat(subColor{isub},sum(sessDistVec(:)>0),1),scatterFill);
                    scatter1.MarkerFaceAlpha = surfaceAlpha;
                    scatter1.MarkerEdgeAlpha = surfaceAlpha;
                end
                hold on
                for idist=1:minSessions-1
                    subMeanDistData(isub,idist) = mean(subData(sessDistVec(:)==idist));
                end
                if singleSubject
                    plot(subMeanDistData(isub,:),'linewidth', linewidthNarrow,'color',subColor{isub}); hold all
                end
                pSubPermDistCorr(isub,i) = sum(subDistCorrPerm(i,isub,:)<=subDistCorr(isub,i))/nperms;
            end
            distCorr(i) = mean(subDistCorr(:,i),1);%mean over subjects
            distCorrPerm(i,:) = mean(subDistCorrPerm(i,:,:),2);%mean over subjects
            pPermDistCorr(i) = sum(distCorrPerm(i,:)<=distCorr(i))/nperms;
            plot(mean(subMeanDistData,1),'linewidth', linewidthWide,'color','k'); hold all
            if ~singleSubject
                dsErrorsurface(1:size(subMeanDistData,2),mean(subMeanDistData,1),std(subMeanDistData)/sqrt(nsubjects),errorbarColor,surfaceAlpha);
            end
    end
    
    if ~saveFigs
        title(['r=' num2str(distCorr(i),'%4.2f') ' p=' num2str(pPermDistCorr(i),'%4.3f')]);
    end
    
    axis square
    xlabel('\Delta session');
    ylabel(ylabelStr);
    xlim([1 size(distData,2)]);
    xticks([1 size(distData,2)]);
    xticklabels([1 size(distData,2)]);
 
    %permutation histogram
    subplot(rows,cols,3+(i-1)*2*cols);
    h=histogram(distCorrPerm(i,:),histBins,'faceColor',faceColor,'edgeColor',edgeColor); hold all
    h.FaceColor = 0.7*[1 1 1];
    h.Normalization = 'probability';
    axis square
    
    xlim([-0.2 0.2]);
    if abs(distCorr(i))>0.2
        xlim([-0.3 0.3]);
    end
    ylmt = get(gca,'ylim');
    plot([distCorr(i) distCorr(i)], [ylmt(1) ylmt(2)],'k','linewidth',linewidthWide);
    xlabel('correlation(r)');
    ylabel('permutations prob.');
    
    if i==1%only when using R2, not Pearson's r
        %generalization to adjacent sessions
        subplot(rows,cols,cols+2); hold all
        idist=1;
        switch coreyVersion
            case 0
                meanDistSess = distSessData{iroi,idist}(subjects,1:minSessions-idist);
                if singleSubject
                    for isub=1:length(subjects)
                        plot(meanDistSess(isub,:),'color',subColor{subjects(isub)},'linewidth',linewidthNarrow); hold all
                    end
                else
                    dsErrorsurface(1:size(meanDistSess,2),mean(meanDistSess),std(meanDistSess)/sqrt(nsubjects),[0.6 0.6 0.6] ,surfaceAlpha);
                end
                plot(mean(meanDistSess),'linewidth',linewidthWide,'color','k'); hold on
                %get p-value for adjacent generalization
                meanDistSessPerm = squeeze(sum(distSessDataPerm));%sum over subjects
                [distSessCorr(i), pDistSessCorr(i)] = corr([1:minSessions-1]',mean(meanDistSess)');
                for iperm=1:nperms
                    [distSessCorrPerm(i,iperm), pDistSessCorrPerm(i,iperm)] = corr([1:minSessions-1]',meanDistSessPerm(iperm,:)');
                end
                pPermDistSessCorr(i) = sum(distSessCorrPerm(i,:)<=distSessCorr(i))/nperms;
            case 1
                %compute correlation using all elements of the goodness-of-fit matrix:
                dataDistMat = squeeze(mean(similarityData,1));%mean over subjects
                dataDistSess = dataDistMat(distMatrix(1:minSessions,1:minSessions)==idist);
                [trainSess testSess] = ind2sub(size(distMatrix(1:minSessions,1:minSessions)),find(distMatrix(1:minSessions,1:minSessions)==idist));%make sure this is correct
                [distSessCorr(i), pDistSessCorr(i)] = corr(trainSess,dataDistSess);
                
                permDistMat = squeeze(mean(similarityPerm,1)); %mean over subjects, for each permutation
                for iperm=1:nperms
                    tempPermMat = squeeze(permDistMat(iperm,1:minSessions,1:minSessions));
                    permDistVec = tempPermMat(distMatrix(1:minSessions,1:minSessions)==idist);
                    [distSessCorrPerm(i,iperm), pDistCorrPerm(i,iperm)] = corr(trainSess,permDistVec);
                end
                pPermDistSessCorr(i) = sum(distSessCorrPerm(i,:)<=distSessCorr(i))/nperms;
                %plot
                if addScatter
                    scatter1=scatter(trainSess,dataDistSess,scatterSize,[0 0 0],scatterFill);
                    scatter1.MarkerFaceAlpha = 1;
                    scatter1.MarkerEdgeAlpha = 1;
                end
                for ifirst=1:minSessions-1
                    meanDistSessData(ifirst) = mean(dataDistSess(trainSess==ifirst));
                end
                plot(meanDistSessData,'linewidth', linewidthWide,'color','k'); hold all
                
            case 2
                %compute correlation per subject and then average
                for isub=1:length(subjects)
                    subData = similarityData(isub,1:minSessions,1:minSessions);
                    dataDistSess = subData(distMatrix(1:minSessions,1:minSessions)==idist);
                    [trainSess testSess] = ind2sub(size(distMatrix(1:minSessions,1:minSessions)),find(distMatrix(1:minSessions,1:minSessions)==idist));%make sure this is correct
                    [subDistSessCorr(isub,i) pSubDistSessCorr(isub,i)] = corr(trainSess,dataDistSess);
                    
                    permDistMat = squeeze(similarityPerm(isub,:,:,:));
                    for iperm=1:nperms
                        tempPermMat = squeeze(permDistMat(iperm,1:minSessions,1:minSessions));
                        permDistVec = tempPermMat(distMatrix(1:minSessions,1:minSessions)==idist);
                        [subDistSessCorrPerm(i,isub,iperm), pSubDistSessCorrPerm(i,isub,iperm)] = corr(trainSess,permDistVec);
                    end
                    %plot
                    if addScatter
                        scatter1=scatter(trainSess,dataDistSess,scatterSize,repmat(subColor{isub},length(trainSess),1),scatterFill);
                        scatter1.MarkerFaceAlpha = 1;
                        scatter1.MarkerEdgeAlpha = 1;
                    end
                    for ifirst=1:minSessions-1
                        subMeanDistSessData(isub,ifirst) = mean(dataDistSess(trainSess(:)==ifirst));
                    end
                    plot(subMeanDistSessData(isub,:),'linewidth', linewidthNarrow,'color',subColor{isub}); hold all
                end
                distSessCorr(i) = mean(subDistSessCorr(:,i),1);%mean over subjects
                distSessCorrPerm(i,:) = mean(subDistSessCorrPerm(i,:,:),2);%mean over subjects
                pPermDistSessCorr(i) = sum(distSessCorrPerm(i,:)<=distSessCorr(i))/nperms;
                plot(mean(subMeanDistSessData,1),'linewidth', linewidthWide,'color','k'); hold all
                for isub=1:length(subjects)
                    pSubPermDistSessCorr(isub,i) = sum(subDistSessCorrPerm(i,isub,:)<=subDistSessCorr(isub,i))/nperms;
                end
        end
        
        axis square
        xlabel('session');
        ylabel(ylabelStr);
        xlim([1 minSessions-idist]);
        xticks([1 minSessions-idist]);
        xticklabels([1 minSessions-idist]);
        if ~saveFigs
            title(['r=' num2str(distSessCorr(i),'%4.2f') ' p=' num2str(pPermDistSessCorr(i),'%4.3f')]);
        end
        
        %adjacent permutation histogram
        subplot(rows,cols,cols+3)
        h=histogram(distSessCorrPerm(i,:),histBins,'faceColor',faceColor,'edgeColor',edgeColor); hold all
        h.FaceColor = 0.7*[1 1 1];
        h.Normalization = 'probability';
        axis square
        ylmt = get(gca,'ylim');
        plot([distSessCorr(i) distSessCorr(i)], [ylmt(1) ylmt(2)],'k','linewidth',linewidthWide);
        xlim([-0.2 0.2]);
        xlabel('correlation(r)');
        ylabel('permutations prob.');
        
        %schematic of adjacent R2
        subplot(rows,cols,cols+1)
        img=imagesc(abs(sessDiff));
        img.AlphaData = abs(sessDiff)==1;
        AlphaDataMapping = 'scaled';

        xlabel('test session');
        ylabel('train session');
                xticks([1 minSessions]);
        yticks([1 minSessions]);
        if ~saveFigs
            title('distance=1');
        end
    end
end

for isubplot=1:rows*cols
    s=subplot(rows,cols,isubplot);
    colormap(s, colormapName);
    axis square
    set(gca, 'box', 'on', 'Visible', 'on');
end
%%
%p-values fpr single subjects
pSubPermDistCorr


for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    set(gca, 'box', 'on', 'Visible', 'on');
end

set(gcf,'position',[250 400 160*cols 600]);
if addColorbars
    set(gcf,'position',[250 300 160*cols 800]);
end
for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    legend off
end
if saveFigs
    if ~addColorbars
        savepdf(f,fullfile(figsFolder,['fig1' fixedFirstStr zscoreStr r2threshStr '.pdf']));
    else
        savepdf(f,fullfile(figsFolder,['fig1' fixedFirstStr zscoreStr r2threshStr '_colorbar.pdf']));
    end
end
