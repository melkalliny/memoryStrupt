function [null, meanCorrs] = plotAggregateCorrelations(plotType,corrType,responseOrProbe,correlations_testGroup,subjects, plotIndividual,plotGroup)
null=[];
meanCorrs = [];

writeOut_folders = {'probe_locked', 'response_locked'};
subjectNames = {'NIH043', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055'};
correlations_testGroup = correlations_testGroup(:,:,responseOrProbe);


timeseries_aggregate = cell(6,5,3,2);
shuffle_aggregate = cell(6,5,3,2);
locsToDo = 1:5;

rawOrNormalized = 'normalized';
locs = {'MTL', 'ATL', 'PTL', 'Ventral Frontal', 'Dorsal Frontal', 'PPC'};


for pt = [subjects]
    for locationToDo = locsToDo
        if plotIndividual==1
            gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
            set(gcf,'Color','w')
        end
        
        maxYLim3 = [NaN, NaN];
        for j=plotType
            if responseOrProbe==2
                        subplotLocs = [1 2 4 3;5 6 8 7];

            else
            subplotLocs = [1 2 5 3 4;6 7 10 8 9];    
            end
            for freqsToDo=1:2
                params = [locationToDo freqsToDo j]; %locations, freqs, corrType
                % get reaction time, to add to plots
                responseTime = correlations_testGroup{params(1),pt}{1,5};
                responseTime = responseTime / 1000; responseTime = responseTime*20;
                timeseries = correlations_testGroup{params(1),pt}{1,params(3)}(params(2),:,1);
                timeseries_aggregate{pt,locationToDo,j,freqsToDo} = timeseries;
                
                % errorbars for real data
                shuffleMean = correlations_testGroup{params(1),pt}{1,params(3)}(params(2),:,2);
                shuffleStd = correlations_testGroup{params(1),pt}{1,params(3)}(params(2),:,3);
                
                shuffle_aggregate{pt,locationToDo,j,freqsToDo} = shuffleMean;
                
                if plotIndividual==0
                    continue
                end
                
                 if responseOrProbe==2
                    windowsToPlot = 30:80;
                    xValues = 1:51;
                else
                    windowsToPlot = 1:51;
                    xValues = 1:51;
                 end

                
                if strcmpi(rawOrNormalized,'normalized')
                    timeseries = timeseries / max(abs(timeseries));
                end
                
                if responseOrProbe==2
                subplot(2,4,subplotLocs(freqsToDo,j))
                else
                subplot(2,5,subplotLocs(freqsToDo,j))   
                end
                plot(xValues,timeseries(windowsToPlot),'k')
                hold on
                
                tempHigh = repmat(nanmean(shuffleMean(windowsToPlot)),1,length(windowsToPlot)) +(1*nanmean(shuffleStd(windowsToPlot)));
                tempLow = repmat(nanmean(shuffleMean(windowsToPlot)),1,length(windowsToPlot)) - (1*nanmean(shuffleStd(windowsToPlot)));

                
                if strcmpi(rawOrNormalized,'normalized')
                    tempHigh = tempHigh / max(abs(timeseries));
                    tempLow = tempLow / max(abs(timeseries));
                end
                
                
                plot(xValues,tempLow,'r')
                hold on
                plot(xValues,tempHigh,'r')
                
                ylimits = ylim;         

                if responseOrProbe==2
                    hold on
                    plot([40 40],[-50 50],'b')
                    xticks([0 20 40])
                    xticklabels([-2 -1 0])
                    if strcmpi(rawOrNormalized,'normalized')
                    ylabel('correlation, normalized')     
                    else
                    ylabel('spearmans corr')
                    end
                    xlabel('time (s)')
                else
                    xticks([0 10 30 50 70])
                    xticklabels([-0.5 0 1 2 3])
                    if strcmpi(rawOrNormalized,'normalized')
                    ylabel('correlation, normalized')     
                    else
                    ylabel('spearmans corr')
                    end
                    xlabel('time (s)')
                end
               
                                    
                if strcmpi(rawOrNormalized,'normalized')
                    ylim([-1.2 1.2])
                else
                    ylim([ylimits(1) ylimits(2)])
                end
                
                                
                maxAbs = max(abs(ylim));
                if j==4; maxYLim3(freqsToDo) = maxAbs; end
                if j==5; maxAbs = maxYLim3(freqsToDo); end
                ylim([-maxAbs maxAbs])
                
                title(sprintf('%s, %d subjs',corrType{j}));
                set(gca,'FontSize',12)
                set(gca, 'FontName', 'Helvetica');
                %set(gca,'FontSize',16)
            end
            
        end
        if plotIndividual==0; continue; end
        suptitle(sprintf('Subject %s, %s',subjectNames{pt},locs{locationToDo}));
        
        if length(find(~isnan(timeseries)))==0 %#ok<*ISMT>
            continue
        else
            writeFigsHere = '/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/figures/ContentVsContext_4plots_individual_updated_without055/';
            if ~exist(sprintf('%s/%s/%s/%s',writeFigsHere,rawOrNormalized,writeOut_folders{responseOrProbe},subjectNames{pt})); mkdir(sprintf('%s/%s/%s/%s',writeFigsHere,rawOrNormalized,writeOut_folders{responseOrProbe},subjectNames{pt})); end
            fig2pngSimple(gcf,sprintf('%s%s/%s/%s/%s.png',writeFigsHere,rawOrNormalized,writeOut_folders{responseOrProbe},subjectNames{pt},locs{locationToDo}));
        end
        
    end
end

meanCorrs = NaN(6,6,2,91);

if plotGroup==1
    for locationToDo = locsToDo
        gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
        set(gcf,'Color','w')
        
        if responseOrProbe==2
            subplotLocs = [1 2 4 3;5 6 8 7];
        else
            subplotLocs = [1 2 5 3 4;6 7 10 8 9];
        end
        
        for j=plotType
            maxYlim = [];
            for freqsToDo=[1,2]
                
                if responseOrProbe==2
                    subplot(2,4,subplotLocs(freqsToDo,j))
                else
                    subplot(2,5,subplotLocs(freqsToDo,j))
                end
                
                all_timeseries = timeseries_aggregate(:,locationToDo,j,freqsToDo);
                all_timeseries = cell2mat(all_timeseries);
                mean_timeseries = all_timeseries;
                
                for subject=1:size(mean_timeseries,1)
                    mean_timeseries(subject,:) = mean_timeseries(subject,:) / max(mean_timeseries(subject,:));
                    
                    if j==4
                     meanCorrs(subject,locationToDo,freqsToDo,:) = mean_timeseries(subject,:);
                    end
                                          
                  if strcmpi(rawOrNormalized,'normalized')
                      all_timeseries(subject,:) = all_timeseries(subject,:) / max(abs(all_timeseries(subject,:)));
                  end
                end
                
                ste_timeseries = nanstd(mean_timeseries) / sqrt(size(mean_timeseries,1));
 
                
                
                if responseOrProbe==2
                    windowsToPlot = 30:80;
                    xValues = 1:51;
                else
                    windowsToPlot = 1:71;
                    xValues = 1:71;
                end
                
                easy_shade_bar_v2(all_timeseries(:,windowsToPlot),'b')
                
                tempMax = max(abs(all_timeseries(:)));
                
                % errorbars for shuffle data
                shuffle_timeseries = shuffle_aggregate(:,locationToDo,j,freqsToDo);
                shuffle_timeseries = cell2mat(shuffle_timeseries);
                hold on
                easy_shade_bar_v2(shuffle_timeseries(:,windowsToPlot),'r')
                
                hold on
                for subject=1:size(all_timeseries,1)
                    tp = plot(xValues,all_timeseries(subject,windowsToPlot),'k');
                    tp.Color(4) = 0.2;
                end
               
                % run t-test on pre-stim vs every window, and put *** above
                % significant time points
                
                if strcmpi(rawOrNormalized,'normalized')
                    ylim([-1.2 1.2])
                else
                    ylim([-tempMax tempMax])
                end
                
                maxAbs = max(abs(ylim));
                if freqsToDo==1; maxYlim(freqsToDo) = maxAbs; end
                if freqsToDo==2; ylim([-maxYlim maxYlim]); end
                
                if responseOrProbe==2
                    ylimits = ylim;
                    hold on
                    plot([40 40],[ylimits(1) ylimits(2)],'b')
                    xticks([0 20 40])
                    xticklabels([-2 -1 0])
                    if strcmpi(rawOrNormalized,'normalized')
                    ylabel('correlation, normalized')     
                    else
                    ylabel('spearmans corr')
                    end
                    xlabel('time (s)')
                else
                    xticks([0 10 30 50 70])
                    xticklabels([-0.5 0 1 2 3])
                    if strcmpi(rawOrNormalized,'normalized')
                    ylabel('correlation, normalized')     
                    else
                    ylabel('spearmans corr')
                    end
                    xlabel('time (s)')
                end
                
                if responseOrProbe==1 && j==1
                   for time=1:length(xValues)
                      temp_test = all_timeseries(:,windowsToPlot(time));
                      temp_baseline = nanmean(all_timeseries(:,1:10),2);
                      [a,b,c,d] = ttest2(temp_test,temp_baseline);
                      if d.tstat > 2
                          hold on
                          ylimits = ylim;
                          plot(windowsToPlot(time),ylimits(2)-0.05*ylimits(2),'o','Color','k')
                      end
                   end
                    
                end

                
                title(sprintf('%s, %d subjs',corrType{j}, length(find(~isnan(all_timeseries(:,1))))));
                set(gca,'FontSize',12)
                set(gca, 'FontName', 'Helvetica');
                
                
            end
        end
        suptitle(sprintf('%s',locs{locationToDo}))
        
        if length(find(~isnan(shuffle_timeseries)))==0 %#ok<*ISMT>
            continue
        else
            writeFigsHere = '/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/figures/ContentVsContext_4plots_updated_without055/';
            if ~exist(sprintf('%s/%s/%s/',writeFigsHere,rawOrNormalized,writeOut_folders{responseOrProbe})); mkdir(sprintf('%s/%s/%s/',writeFigsHere,rawOrNormalized,writeOut_folders{responseOrProbe})); end
            fig2pngSimple(gcf,sprintf('%s%s/%s/%s.png',writeFigsHere,rawOrNormalized,writeOut_folders{responseOrProbe},locs{locationToDo}));
        end
        
        
    end
   % suptitle(sprintf('All Subjects, %s',locationStrings{1,params(1)}));
    
    meanCorrs = meanCorrs(:,:,:,windowsToPlot);
    
end





end

