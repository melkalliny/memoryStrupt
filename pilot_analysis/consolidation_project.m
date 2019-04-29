


%% analysis to-do notes
% window-size exploration
% create new model, for correct-vs-incorrect, with incorrect being double-incorrects
% run correct-vs-incorrect ERP analysis
% run correct-vs-incorrect power analysis



%% does z-scored power change over sessions
consolidation_computePower
consolidation_preprocessPower
consolidation_visualizePower



%% same analyses as above, but for ripples
consolidation_computeRipples
consolidation_preprocessRipples
consolidation_visualizeRipples



%% data-driven approach, looking at correct-vs-incorr and use that
consolidation_preprocessPowerAllFrequencies
[pointsOfInterest] = consolidation_visualizePower_allFrequencies;
consolidation_calculateSimilarity(pointsOfInterest);



%% consistency vs reaction time, group level analysis
freqsOfInterest = 3:17; %60-70, indexing the variable 'freqs'
locationStrings = {'Hippocampus', 'Entorhinal', 'Lateral temporal', 'Frontal', 'Parahippocampal'};
subjectNames = {'NIH049', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055'};
subjectsToRun = [1:6];
[correlations_testGroup] = consolidation_calculateSimilarity_retrieval(subjectNames, subjectsToRun, freqsOfInterest); 
subjectNames = {'NIH029', 'NIH032', 'NIH036', 'NIH041', 'NIH042', 'NIH046', 'NIH047'};
subjectsToRun = 1:7;
[correlations_controlGroup] = consolidation_calculateSimilarity_retrieval(subjectNames, subjectsToRun, freqsOfInterest);
region = 1;
plotGroupConsistencyCorrelations(correlations_testGroup, [], locationStrings, region)



%% all model-based correlation analyses
pool = parpool('local',8);
consolidation_computePower
consolidation_preprocessLowAndHigh
subjectNames = {'NIH043', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055'};

correlations_results = cell(1,12);
subjectsToRun = [1:6];
for doParam = [1:6]
    strings = {'PROBE_START','RESPONSE'};
    params = [1 2 3 4 5 6 1 2 3 4 5 6; 2 2 2 2 2 2 1 1 1 1 1 1];
    stringEvent = strings{params(2,doParam)};
    [correlations_results{1,doParam}] = consolidation_calculateSimilarity_LowAndHigh(subjectNames, subjectsToRun(params(1,doParam)),stringEvent);
end
save('/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/correlations_response_newRT.mat','correlations_results')


parfor doParam = [7:12]
    strings = {'PROBE_START','RESPONSE'};
    params = [1 2 3 4 5 6 1 2 3 4 5 6; 2 2 2 2 2 2 1 1 1 1 1 1];
    stringEvent = strings{params(2,doParam)};
    [correlations_results{1,doParam}] = consolidation_calculateSimilarity_LowAndHigh(subjectNames, subjectsToRun(params(1,doParam)),stringEvent);
end
save('/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/correlations_full_all.mat','correlations_results')


%% plot results of above analyses
plotType = 1:5;
corrType = {'Context', 'Content', 'RT', 'Content-Correct', 'Content-Incorrect'};
responseOrProbe = 1;
plotIndividual = 0; plotGroup = 1;
subjects = [1:5];
correlations_aggregate = cell(6,6,2);
for j=1:6
    for type=1:2
    correlations_aggregate(:,subjectsToRun(j),type) = correlations_results{1,j}(:,subjectsToRun(j));
    end
end
[tROIS, meanCorrs] = plotAggregateCorrelations(plotType,corrType,responseOrProbe,correlations_aggregate,subjects,plotIndividual,plotGroup);


readme = 'first dimension=subjects,second dimension=region, third dimension=lowOrHighFreq, fourth dimension=windows';
save('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/regression_inputs_meanCorrs.mat','meanCorrs', 'readme')



%% export data for the regression based on tROIs
regions = {'MTL', 'ATL', 'PTL', 'V-Frontal', 'D-Frontal', 'PPC'};
regression_inputs = cell(1,1);
subjectsToRun = [1:6];
params(1,1) = 2; params(1,2) = 4; params(1,3) = 1; params(1,4) = 2; % ventral-frontal
params(2,1) = 2; params(2,2) = 4; params(2,3) = 2; params(2,4) = 2; % ventral-frontal
params(3,1) = 2; params(3,2) = 5; params(3,3) = 1; params(3,4) = 2; % dorsal-frontal
params(4,1) = 2; params(4,2) = 5; params(4,3) = 2; params(4,4) = 2; % dorsal-frontal
params(5,1) = 2; params(5,2) = 1; params(5,3) = 2; params(5,4) = 2; % MTL, response
params(6,1) = 2; params(6,2) = 2; params(6,3) = 1; params(6,4) = 2; % ATL, response
params(7,1) = 2; params(7,2) = 2; params(7,3) = 2; params(7,4) = 2; % ATL, response
params(8,1) = 2; params(8,2) = 3; params(8,3) = 2; params(8,4) = 2; % PTL, response
% probe or response locked, which region, low or high freq, how many s
% before response to use
for set = 1:size(params,1)
[regression_inputs{1,set}] = consolidation_extractSimilarity_LowAndHigh(subjectNames, [1:6], params(set,:));
end
save('/Volumes/Shares/FRNU/dataWorking/consolidationProject/intermediates/regression_inputs_updated_v4.mat','regression_inputs')






%% contrast reaction time changes in groups of subjects 
testGroup = {'NIH043', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055'};
%testGroup = {'NIH055'};
controlGroup = {'NIH029', 'NIH032', 'NIH036', 'NIH041', 'NIH042', 'NIH046', 'NIH047'};
groupLevel_repetitionVsRT(testGroup,controlGroup)





%% explore effect of window-size choice on analysis

consolidation_preprocessLowAndHigh_allTime

windowLength = [100 250 500 1000 2000];
windowShift = [10 25 50 100 200];
for windowOptions = 1:3
    subjectNames = {'NIH043', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055'};
    params = [1 2 3 4 5 6 1 2 3 4 5 6; 2 2 2 2 2 2 1 1 1 1 1 1];
    
    correlations_results = cell(1,12);
    subjectsToRun = [1:6];
    for doParam = [4]
        strings = {'PROBE_START','RESPONSE'};
        params = [1 2 3 4 5 6 1 2 3 4 5 6; 2 2 2 2 2 2 1 1 1 1 1 1];
        stringEvent = strings{params(2,doParam)};
        
   
        [correlations_results{1,doParam}] = consolidation_calculateSimilarity_LowAndHigh(subjectNames, subjectsToRun(params(1,doParam)),stringEvent, windowOptions);
        
    end
    %save(sprintf('/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/correlations_windowSizes/windowSize%d.mat',windowOptions),'correlations_results')
    
end




consolidation_createSwarmScripts





