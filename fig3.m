% fig3.m
%
% associated with the following publication: Roth, ZN, and Merriam, EP (2023).
% Representations in human primary visual cortex drift over time
% DOI:
%
%   usage: fig3()
%   by: zvi roth
%   date: 3/10/2022
%   purpose:  Plots for Fig 3
%   uses files created by: savePermPop_expand.m
%   creates pdf files  

close all
clear all
tic

saveFigs=0;
addColorbars = 0;

addScatter = 0;
nperms = 1000;
toZscore = 0;%0, 1, 2==zeroMean, 3==equalStd
fixedFirst = 0;
singleSubject=1;

figRoi = 1;
subjects = 1:8;%[5:8];
nsubjects = length(subjects);

saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand');
figsFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand_figs');
if ~isfolder(saveFolder)
    saveFolder = ['/misc/data18/rothzn/nsd/repDrift_expand/'];
    figsFolder = ['/misc/data18/rothzn/nsd/repDrift_expand_figs/'];
end

histBins = 20;
scatterSize=5;
scatterFill = 'filled';

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

r2thresh = 0;
r2threshStr = '';
if r2thresh>0
    r2threshStr = ['r2thresh' num2str(r2thresh,'%4.2f')];
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
fixedFirstStr='';
if fixedFirst
    fixedFirstStr = '_fixedFirst_';
end
load(fullfile(saveFolder, ['permPop' fixedFirstStr num2str(nperms) zscoreStr r2threshStr  '.mat']),...
    'rois','toNormalize','toZscore','r2thresh','nrois','version', ...
    'permOrders', 'subSessions', 'subjects', 'minSessions','distMatrix',...
    'betweenSessCorr','betweenSessCorrOri',...
    'avgImgCorrMat','avgImgCorrMatOri',...
    'betweenSessImg','betweenSessImgOri','betweenSessImgPerm','betweenSessImgOriPerm',...
    'betweenSessDist','betweenSessDistOri','betweenSessDistPerm','betweenSessDistOriPerm',...
    'numGoodVox',...
    'imgCorrMatPerm','imgCorrMatOriPerm','betweenSessCorrPerm','betweenSessCorrOriPerm');
%%
sessDistVec = distMatrix(:);

%% distance matrix
sess = 1:minSessions;
sessDiff = sess - sess';
%% only chosen ROI - correlation between RDMs
ifig=1;
for coreyVersion=2:2
    figure(ifig); ifig=ifig+1; clf
    iroi=figRoi;
    rows=3;
    cols=2;
    for i=1:cols
        switch i
            case 1
                corrMat = avgImgCorrMatOri;
                distVec = betweenSessImgOri;
                distVecPerm = betweenSessImgOriPerm;
                titleStr = 'pop resp';
                corrMatPerm = imgCorrMatOriPerm;
            case 2
                corrMat = betweenSessCorrOri;
                distVec = betweenSessDistOri;
                distVecPerm = betweenSessDistOriPerm;
                titleStr = 'RDMs';
                corrMatPerm = betweenSessCorrOriPerm;
                
        end
        subplot(rows,cols,i)
        temp = squeeze(nanmean(corrMat(iroi,:,1:minSessions,1:minSessions),2));
        img=imagesc(temp,'AlphaData',abs(sessDiff)>0); axis square%mean across subjects
        caxis([min(temp(abs(sessDiff)>0)), max(temp(abs(sessDiff)>0))]);
    
        
        if addColorbars
            colorbar('southoutside')
            xticks([]);
            yticks([]);
        else 
            xlabel('session'); ylabel('session');
            xticks([1 minSessions]);
            yticks([1 minSessions]);
        end
        if ~saveFigs
            title(titleStr);
        end
        subplot(rows,cols,i+cols)
        
        
        switch coreyVersion
            case 0
                if singleSubject
                    for isub=1:length(subjects)
                       plot(squeeze(distVec(iroi,isub,:)),'color',graySubColor,'linewidth',linewidthNarrow); hold all
                    end
                else
                    dsErrorsurface(1:size(distVec,3),squeeze(mean(distVec(iroi,:,:),2)),squeeze(std(distVec(iroi,:,:),0,2))/sqrt(nsubjects),errorbarColor,surfaceAlpha); hold on
                end
                meanDist = squeeze(mean(distVec(iroi,:,:),2));
                plot(meanDist,'linewidth',linewidthWide,'color','k'); axis square
                meanDistPerm = squeeze(mean(distVecPerm(iroi,:,:,:),2));
                [distCorr(i), pDistCorr(i)] = corr([1:minSessions-1]',meanDist);
                for iperm=1:nperms
                    [distCorrPerm(i,iperm), pDistCorrPerm(i,iperm)] = corr([1:minSessions-1]',meanDistPerm(iperm,:)');
                end
                pPermDistCorr(i) = 1-sum(distCorrPerm(i,:)>=distCorr(i))/nperms;
                
            case 1
                %compute correlation using all elements of the goodness-of-fit matrix:
                dataDistMat = squeeze(mean(corrMat(iroi,:,1:minSessions,1:minSessions),2));%mean over subjects
                [distCorr(i), pDistCorr(i)] = corr(distMatrix(distMatrix>0),dataDistMat(distMatrix>0));
                permDistMat = squeeze(mean(corrMatPerm(iroi,:,:,:,:),2)); %mean over subjects, for each permutation
                for iperm=1:nperms
                    permDistVec = permDistMat(iperm,1:minSessions,1:minSessions);
                    [distCorrPerm(i,iperm), pDistCorrPerm(i,iperm)] = corr(distMatrix(distMatrix>0),permDistVec(distMatrix>0));
                end
                pPermDistCorr(i) = 1-sum(distCorrPerm(i,:)>=distCorr(i))/nperms;
                %plot
                if addScatter
                    scatter1=scatter(distMatrix(distMatrix>0),dataDistMat(distMatrix>0),scatterSize,[0 0 0],scatterFill); hold all
                    scatter1.MarkerFaceAlpha = 2*surfaceAlpha;
                    scatter1.MarkerEdgeAlpha = 2*surfaceAlpha;
                end
                for idist=1:minSessions-1
                    meanDistData(idist) = mean(dataDistMat(distMatrix(:)==idist));
                end
                plot(meanDistData,'linewidth', linewidthWide,'color','k'); 
                
            case 2
                %compute correlation per subject and then average
                for isub=1:length(subjects)
                    subData = corrMat(iroi,isub,1:minSessions,1:minSessions);
                    [subDistCorr(isub,i), pSubDistCorr(i)] = corr(distMatrix(distMatrix>0),subData(distMatrix>0));
                    permDistMat = squeeze(corrMatPerm(iroi,isub,:,:,:)); %this subjects, all permutations
                    for iperm=1:nperms
                        permDistVec = permDistMat(iperm,1:minSessions,1:minSessions);
                        [subDistCorrPerm(i,isub,iperm), pSubDistCorrPerm(i,isub,iperm)] = corr(distMatrix(distMatrix>0),permDistVec(distMatrix>0));
                    end
                    %plot
                    if addScatter
                        scatter1 = scatter(distMatrix(distMatrix>0),subData(distMatrix>0),scatterSize,repmat(subColor{isub},sum(distMatrix(:)>0),1),scatterFill);
                        scatter1.MarkerFaceAlpha = surfaceAlpha;
                        scatter1.MarkerEdgeAlpha = surfaceAlpha;
                    end
                    hold on
                    for idist=1:minSessions-1
                        subMeanDistData(isub,idist) = mean(subData(distMatrix(:)==idist));
                    end
                    if singleSubject
                        plot(subMeanDistData(isub,:),'linewidth', linewidthNarrow,'color',subColor{isub}); hold all
                    end
                    pSubPermDistCorr(isub,i) = 1-sum(subDistCorrPerm(i,isub,:)>=subDistCorr(isub,i))/nperms;
                end
                distCorr(i) = mean(subDistCorr(:,i),1);%mean over subjects
                distCorrPerm(i,:) = mean(subDistCorrPerm(i,:,:),2);%mean over subjects
                pPermDistCorr(i) = 1-sum(distCorrPerm(i,:)>=distCorr(i))/nperms;
                plot(mean(subMeanDistData,1),'linewidth', linewidthWide,'color','k'); hold all
                if ~singleSubject
                    cla
                    temp = subMeanDistData;
                    tempDiff = temp - temp(:,1);
                    tempPerc = 100*tempDiff./temp(:,1);
                    plot(median(tempPerc,1),'linewidth', linewidthWide,'color','k'); hold all
                    dsErrorsurface(1:size(tempPerc,2),median(tempPerc,1),std(tempPerc)/sqrt(nsubjects),errorbarColor,surfaceAlpha);

                    plot(1:29,zeros(1,29),'k');
                    ylim([-7 3]);
                end
        end
        if ~saveFigs
            title(['r=' num2str(distCorr(i),'%4.2f') ' p=' num2str(pPermDistCorr(i),'%4.3f')]);
        end
        ylabel('correlation (r)');
        if ~singleSubject
            ylabel('median % change correlation');
        end
        xlabel('\Deltasession');
                xticks([1 minSessions-1]);

        subplot(rows,cols,i+2*cols)
        h=histogram(distCorrPerm(i,:),histBins,'faceColor',faceColor,'edgeColor',edgeColor); hold all
        h.FaceColor = 0.7*[1 1 1];
        h.Normalization = 'probability';
        axis square
        ylmt = get(gca,'ylim');
        plot([distCorr(i) distCorr(i)], [ylmt(1) ylmt(2)],'k','linewidth',linewidthWide);
        xlim([-0.2 0.2]);
        xlabel('correlation(r)');
        ylabel('permutations prob.');

    end
    for isubplot=1:rows*cols
        subplot(rows,cols,isubplot)
        set(gca, 'box', 'on', 'Visible', 'on');
        axis square
    end
    set(gcf,'position',[250+coreyVersion*90 200+coreyVersion*20 500 450]);
end
if saveFigs
    if ~addColorbars
        print('-painters','-dpdf',fullfile(figsFolder,['brief3both' zscoreStr  '.pdf']));
    else
        print('-painters','-dpdf',fullfile(figsFolder,['brief3both' zscoreStr  '_colorbar.pdf']));
    end
end

pSubPermDistCorr
toc