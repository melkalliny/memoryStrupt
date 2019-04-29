function [null] = groupLevel_repetitionVsRT(testGroup,controlGroup)
null=[];

groups = cell(1,2);
groups{1,1} = testGroup;
groups{1,2} = controlGroup;

RTdiff_testGroup = NaN(3,length(testGroup));
RTdiff_controlGroup = NaN(3,length(testGroup));

 meanCorrectBySession = NaN(6,6);
 

for groupIndex=1
    for n=1:length(groups{1,groupIndex})
        if groupIndex==1
            load(sprintf('/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/processedPower_LowAndHigh_monopolar/%s/allPowData.mat',groups{1,groupIndex}{1,n}),'evsOfInterest')
        elseif groupIndex==2
            load(sprintf('/Volumes/Shares/FRNU/data/eeg/%s/behavioral/palRam/events.mat',groups{1,groupIndex}{1,n}));
            [events(find(strcmp({events.type},'STUDY_PAIR'))).type] = deal('STUDY_PAIR_START');
            [events(find(strcmp({events.type},'TEST_PROBE'))).type] = deal('PROBE_START');
            evsPull = events((union(find(strcmp({events.type},'STUDY_PAIR_START')),find(strcmp({events.type},'PROBE_START')))));
            
            % add events centered around reaction time
            evsPull_plusRT = [];
            for ev=1:length(evsPull)
                evsPull_plusRT = cat(1,evsPull_plusRT,evsPull(ev));
                if isnumeric(evsPull(ev).RT) && strcmpi(evsPull(ev).type,'PROBE_START')
                    newElem = evsPull(ev);
                    newElem.eegoffset = evsPull(ev).eegoffset + evsPull(ev).RT;
                    newElem.type = 'RESPONSE';
                    evsPull_plusRT = cat(1,evsPull_plusRT,newElem);
                elseif ~isnumeric(evsPull(ev).RT) && strcmpi(evsPull(ev).type,'PROBE_START')
                    fprintf('\nno annotated RT - adding event for RT here\n')
                    newElem = evsPull(ev);
                    newElem.eegoffset = evsPull(ev).eegoffset + 0;
                    evsPull_plusRT = cat(1,evsPull_plusRT,newElem);
                    %keyboard
                end
                
            end
            evsPull = evsPull_plusRT;
            evsOfInterest = evsPull(find(strcmp({evsPull.type},'PROBE_START')));
        end
          sessions = unique([evsOfInterest.session]);
          
          
          evsOfInterest = evsOfInterest(find(strcmp({evsOfInterest.type},'PROBE_START')));
          
          meanRTBySession = NaN(1,length(sessions));
         
          
          for i=1:length(sessions)
              tempEvs = evsOfInterest(find([evsOfInterest.session] == sessions(i)));
              findCorrect = tempEvs(find([tempEvs.correct] == 1));
              correctRTs = [findCorrect.RT];
              correctRTs(find(correctRTs==-999)) = NaN;
              meanRTBySession(i) = nanmean(correctRTs);   
              meanCorrectBySession(n,i) = length(findCorrect) / length(tempEvs);
          end
          
          [rho,p] = corr((1:length(meanRTBySession))',meanRTBySession','Type','Spearman');
          
          
         
        
        if groupIndex==1
            RTdiff_testGroup(1,n) = rho;
        elseif groupIndex==2
            RTdiff_controlGroup(1,n) = rho;
        end
        
        
        
    end
end


tempTest = RTdiff_testGroup(1,:);
[a,b,c,d] = ttest(tempTest);
tempControl = RTdiff_controlGroup(1,:);
[a,b,c,d] = ttest2(tempTest,tempControl);

    gcf = figure('Units', 'Normalized', 'OuterPosition', [0 0 0.8 0.8]);
    set(gcf,'Color','w')
%bar(1:2,[nanmean(tempTest) nanmean(tempControl)])
testEB = nanstd(tempTest) / sqrt(6);
testControl = nanstd(tempControl) / sqrt(7);
hold on
errorbar(1:2,[nanmean(tempTest) nanmean(tempControl)],[testEB testControl])
ylim([-1 0])
xlim([0.5 2.5])
xticks([1 2])
ylabel('Rho value')
set(gca, 'FontName', 'Helvetica');
set(gca, 'FontSize', 14);





end

