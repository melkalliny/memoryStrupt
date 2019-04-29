function [iChanKeep,iEvKeep,strOut] = cleanEEGevents(eegWave, weights)
%
%  Pass in a 3D matrix of EEG data... this function will find channels and events
%    with unreasonable variance and kurtosis and trim them out
%
%  input:   eegWave (chan x events x time) -- mono or bipolar. could also be power at one band
%           FIG_NUM   -- figure to plot cleaning result on
%           FIG_TITLE -- text to be added to result figure (e.g., 'NIH030')
%           weights   -- [w_t,w_c]; %default- jw[0.5 0.5], RAM [2.3 2.3]
%
%
%  outputs: index of channels to keep
%           index of events to keep
%
%

fprintf('\nCleaning EEG (channels x trials): ');

% concatenate variance and kurtosis for both data structures
% krtall = [chanVSev(wKUR).data];                                %- must be (channels x trials)
% varall = [chanVSev(wVAR).data];          varStr='Var(EEG)';    F_OFFSET=10;%- var
% %varall = [chanVSev(pVAR).data(:,:,1)];   varStr='Var(LoFpow)'; F_OFFSET=20;%- var
% %varall = [chanVSev(pVAR).data(:,:,end)]; varStr='Var(HiFpow)'; F_OFFSET=30;%- var
% if size(varall,1)~= numChannels | size(varall,2)~=length(eventTrigger), fprintf('\nUh oh. check dimensions of krtall'); keyboard; end


varall    = var(eegWave,0,3);       varStr='Var(EEG)';  
krtall    = kurtosis(eegWave,1,3);           
   
iCH = 1; iTR = 2;  %-channels should be first dimension, trials should be second. If not, while loop below wont work correctly
numChan   = size(eegWave,iCH);
numEvents = size(eegWave,iTR);


% set weight for determining outliers
% default for boxplot is 1.5
% adjust upward to be more lenient, downward to be more aggressive
% the two variables will balance each other; in general, adjusting one will
% cause the other to effectively compensate in the other direction
if nargin<2,
    w_t = 0.5;
    w_c = 0.5;
else
    w_t = weights(1);  % (2.3 RAM) weight for trials    -- 0.5 trims upper 90% of a norm distribution (~1.3 SD); 1.0 trims upper 2.5% (~2 SD); 1.8 trims upper 0.1% (~3 SD);  2.3 trims upper <.001% (~3.7 SD if normal!!!)
    w_c = weights(2);  % (2.3 RAM) weight for channels
end

%- weights are analagous to SD's for a normal distribution, where each 0.5 weight adds another SD to the threshold (which has a minimum value of 1SD, so wt=0.5 --> 2SD threhsold)
%- the non-parametric approximation of SD is the range of the 25th-to-75th percentile, that is, the range of the middle 50% of the data.  
%%%-  run the following code to see how this works with a normal distribution -%%%
% varall=randn(1,100000); %- normal distribution, SD=1
% figure; cdfplot(varall);
% wList=[0:.1:3]; dOut=[];
% for w_t=wList,
%     trlthresh1 = quantile(max(varall,[],iCH),0.75) + w_t*(quantile(max(varall,[],iCH),0.75) - quantile(max(varall,[],iCH),0.25)); %-max across channels to find trials  --- 75th quartile + w_t*(25-to-75 quartile width) = (~1SD) + (~1SD per 0.5 w_t) if normal
%     pcntCut = 100.0*sum(varall>trlthresh1)/length(varall);
%     dOut=[dOut; w_t trlthresh1 pcntCut];
% end
% disp('Weight  thresh (numSD)  pcntCut')
% disp(dOut)



% find outlying trials using max variance for each trial (across channels)
trlthresh1 = quantile(max(varall,[],iCH),0.75) + w_t*(quantile(max(varall,[],iCH),0.75) - quantile(max(varall,[],iCH),0.25)); %-max across channels to find trials  --- 75th quartile + w_t*(25-to-75 quartile width) = (~1SD) + (~1SD per 0.5 w_t) if normal
%trlthresh2 = quantile(max(varall,[],iCH),0.25) - w_t*(quantile(max(varall,[],iCH),0.75) - quantile(max(varall,[],iCH),0.25)); %-max across channels to find trials  --- 25th quartile - w_t*(25-to-75 quartile width) = (~1SD) + (~1SD per 0.5 w_t) if normal
%trlselvardum = find(max(varall,[],iCH)>trlthresh1 || min(varall,[],iCH)<trlthresh2);
trlselvardum = find(max(varall,[],iCH)>trlthresh1);

% find outlying channels using max variance for each channel (across trials)
chnthresh1 = quantile(max(varall,[],iTR),0.75) + w_c*(quantile(max(varall,[],iTR),0.75) - quantile(max(varall,[],iTR),0.25)); %-max across trials to find channels
chnselvardum = find(max(varall,[],iTR)>chnthresh1);

% VS added: to get trials or channels with no data
emptyTrl = find(isnan(mean(krtall,iCH)));
emptyChn = find(isnan(mean(krtall(:,find(~isnan(mean(krtall,iCH)))),iTR)));


trlselvar = [];
chnselvar = [];
varalldum = varall;
while(1)
    while(~isempty(trlselvardum))
        % remove the furthest outlying trial in each iteration
        [~,maxvartrlind] = max(max(varalldum(:,trlselvardum)));
        varalldum(:,trlselvardum(maxvartrlind)) = nan(size(varalldum,1),1);
        chnthreshnew = quantile(max(varalldum,[],2),0.75) + w_c*(quantile(max(varalldum,[],2),0.75) - quantile(max(varalldum,[],2),0.25));
        [~,chnrej] = setdiff(chnselvardum,find(max(varalldum,[],2)>chnthreshnew));
        % if removing this trial changed the number of outlying channels, then
        % set this trial aside for possible permanent exclusion
        if ~isempty(chnrej)
            chnselvardum = chnselvardum(setxor(1:length(chnselvardum),chnrej));
            trlselvar = [trlselvar trlselvardum(maxvartrlind)];
            %disp(['Adding trial #' num2str(trlselvardum(maxvartrlind)) ' to trlselvar'])
        end
        trlselvardum = trlselvardum(setxor(1:length(trlselvardum),maxvartrlind));
    end
    
    % after removing all outlying trials, check to see if any channels are still
    % outlying; if so, set the furthest outlying channel aside for exclusion,
    % then clear the temporary trial exclusion list
    if ~isempty(chnselvardum)
        [~,maxvarchnind] = max(max(varalldum(chnselvardum,:),[],2));
        chnselvar = [chnselvar; chnselvardum(maxvarchnind)];
        %disp(['Adding channel #' num2str(chnselvardum(maxvarchnind)) ' to chnselvar'])
        %disp('Clearing trlselvar')
        trlselvar = [];
    end
    
    % reset the variance measures; keep permanent excluded trials and channels out
    varalldum = varall;
    varalldum(:,trlselvar) = nan(size(varalldum,1),length(trlselvar));
    varalldum(chnselvar,:) = nan(length(chnselvar),size(varalldum,2));
    
    % check again for outlying channels; if there are none, then break
    chnthreshdum = quantile(max(varalldum,[],2),0.75) + w_c*(quantile(max(varalldum,[],2),0.75) - quantile(max(varalldum,[],2),0.25));
    chnselvardum = find(max(varalldum,[],2)>chnthreshdum);
    if isempty(chnselvardum)
        break
    end
    % with outlying channels excluded, repopulate temporary outlying trial list
    trlthreshdum = quantile(max(varalldum),0.75) + w_t*(quantile(max(varalldum),0.75) - quantile(max(varalldum),0.25));
    trlselvardum = find(max(varalldum)>trlthreshdum);
end
trlselvar = sort(trlselvar); % trials to exclude based on variance
chnselvar = sort(chnselvar); % channels to exclude based on variance

krtallnew = krtall;
krtallnew(:,trlselvar) = nan(size(krtallnew,1),length(trlselvar));
krtallnew(chnselvar,:) = nan(length(chnselvar),size(krtallnew,2));

% for kurtosis, only eliminate outlying trials, since we don't
% expect to see entire channels with high kurtosis
w_k = 5; % weight for trials, kurtosis
% find outlying trials using max kurtosis
krtthresh = quantile(max(krtallnew),0.75) + w_k*(quantile(max(krtallnew),0.75) - quantile(max(krtallnew),0.25));
trlselkrt = find(max(krtallnew)>krtthresh);


% trlsel indicates the trials marked for exclusion
trlsel = union(trlselvar,trlselkrt);
% chnsel indicates the channels marked for exclusion
chnsel = chnselvar;

% VS added: eliminate empty recordings
trlsel = sort(union(trlsel, emptyTrl));
chnsel = sort(union(chnsel, emptyChn));


%%- OUTPUT VARIABLES: trials and channels that are cut
strOut    = sprintf(' EVENTS CLEANED: cut %d/%d events, cut %d/%d channels <%s wt t%.1f,c%.1f>',length(trlsel),numEvents,length(chnsel),numChan,varStr,w_t,w_c);
iEvCut    = trlsel;
iChanCut  = chnsel;
iEvKeep   = setdiff([1:numEvents],trlsel);
iChanKeep = setdiff([1:numChan],chnsel);



% optional: plot max variance per trial/channel
PLOT_VAR_REDUCTION=0;
if PLOT_VAR_REDUCTION,
    fSize = 18; %- fontsize
    varallnew = varall;
    varallnew(:,trlsel) = nan(size(varallnew,1),length(trlsel));
    varallnew(chnsel,:) = nan(length(chnsel),size(varallnew,2));
    
    figure(FIG_NUM); clf; set(gcf,'color','w');
    
    subplot(2,4,1); set(gca,'fontsize',fSize);
    scatter(1:size(varall,2),max(varall),'o')
    ax1 = gca;
    ylabel(varStr)
    title(sprintf('%s: %s before rejection',FIG_TITLE,varStr))
    subplot(2,4,6); set(gca,'fontsize',fSize);
    scatter(max(varall,[],2),1:size(varall,1),'o')
    axis ij
    ax2 = gca;
    xlabel(varStr)
    subplot(2,4,5); set(gca,'fontsize',fSize);
    imagesc(varall)
    colormap(jet) %- different matlab versions use different colormaps
    axis ij
    ax3 = gca;
    xlabel('Trial')
    ylabel('Channel')
    xh=xlim; yh=ylim;
    set(ax1,'XTickLabel',[],'TickDir','out','XLim',xh)
    set(ax2,'YTickLabel',[],'TickDir','out','YLim',yh)
    
    subplot(2,4,3); set(gca,'fontsize',fSize);
    scatter(1:size(varallnew,2),max(varallnew),'o')
    ax4 = gca;
    ylabel(varStr)
    title(sprintf('%s after rejection <wt t%.1f,c%.1f>',varStr,w_t,w_c))
    subplot(2,4,8); set(gca,'fontsize',fSize);
    scatter(max(varallnew,[],2),1:size(varallnew,1),'o')
    axis ij
    ax5 = gca;
    xlabel(varStr)
    subplot(2,4,7); set(gca,'fontsize',fSize);
    imagesc(varallnew)
    axis ij
    ax6 = gca;
    xlabel('Trial')
    ylabel('Channel')
    xh=xlim; yh=ylim;
    set(ax4,'XTickLabel',[],'TickDir','out','XLim',xh)
    set(ax5,'YTickLabel',[],'TickDir','out','YLim',yh)
    
    subplot(2,4,4); set(gca,'fontsize',fSize);
    hold on
    for k1 = 1:length(trlsel)
        line([trlsel(k1) trlsel(k1)],yh,'Color','r','LineWidth',1)
    end
    for k1 = 1:length(chnsel)
        line(xh,[chnsel(k1) chnsel(k1)],'Color','r','LineWidth',1)
    end
    xlim(xh);ylim(yh)
    axis ij
    ax7 = gca;
    title([num2str(length(trlsel)) ' trls; ' num2str(length(chnsel)) ' chns'])
    set(gcf,'Position',[68 219 1069 524])
    set(ax1,'Position',[0.05 0.65 0.25 0.27])
    set(ax2,'Position',[0.33 0.1 0.13 0.5])
    set(ax3,'Position',[0.05 0.1 0.25 0.5],'TickDir','out')
    set(ax4,'Position',[0.51 0.65 0.25 0.27])
    set(ax5,'Position',[0.79 0.1 0.13 0.5])
    set(ax6,'Position',[0.51 0.1 0.25 0.5],'TickDir','out')
    set(ax7,'Position',[0.79 0.65 0.13 0.27],'XTickLabel',[],'YTickLabel',[],'TickDir','out')
    clear ax*
end


