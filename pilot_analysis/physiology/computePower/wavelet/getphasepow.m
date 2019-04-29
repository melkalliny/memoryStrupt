function [varargout] = getphasepow(chan,events,DurationMS,OffsetMS,BufferMS,varargin)
%GETPHASEPOW - Calculate wavelet phase and power for a set of events.
%
% Calculate wavelet phase and power as a function of time and frequency for a
% single electrode.
%
% FUNCTION:
%   [phase,pow,kInd] = getphasepow(chan,events,DurationMS,OffsetMS,BufferMS,varargin)
%
% INPUT ARGs:
%   chan = 2;
%   events = events;
%   DurationMS = 2000;
%   OffsetMS = 0;
%   BufferMS = 1000;
%
%   OPTIONAL PARAMS:
%     'freqs'
%     'width'
%     'filtfreq'
%     'filttype'
%     'filtorder'
%     'resampledrate' - resample applied before calculating phase/power
%     'downsample' - decimate applied after calculating power
%     'powonly'
%     'usesingles'
%     'kthresh'- Kurtosis threshold to throw out events.
%
% OUTPUT ARGS:
%   phase- (Events,Freqs,Time)
%   pow- (Events,Freqs,Time)
%   kInd - logical indexes of the events that were not thrown out due to kurtosis
%   (if no phase is demanded then kInd will be the second output argument)

% Changes:
%
% 1/6/06 - PBS - Added decimation of phase.
% 9/15/05 - PBS - Added downsampling via decimate following power calculation.
% 9/15/05 - PBS - Return the logical index of the events not thrown out
%                 with the kurtosis thresh.
% 7/8/05 - MvV - Return the index of the events thrown out with the
%                kurtosis threshold if desired
% 1/18/05 - PBS - Ignore the buffer when applying kurtosis.
%                 Changed round to fix when determine durations.
%                 No longer gets double buffer when buffer is specified.
% 10/30/04 - PBS - Added a kthresh option for filtering events with
%                  high kurtosis.
% 8/25/04 - PBS - Made it an option to create phase and pow
%                 matrixes as singles.
% 3/18/04 - PBS - Switched to gete_ms and added ability to resample
% 11/20 josh 1. changed filtfreq to [] rather than 0. 2. Got rid of the 'single()' function calls that were
% forcing us to return single rather than doubles.  it was pretty
% annoying. per was there a good reason for this function to return
% singles?


% set the defaults
freqs = [];
width = [];
filtfreq =  [];
filttype = 'stop';
filtorder = 1;
resampledrate = [];
dsample = [];
powonly = 0;
usesingles = 0;
keepk = 0;
kthresh = [];
% process the varargs
i = 1;
if( ~isempty(varargin) )
    while i<=length(varargin)
        switch lower(varargin{i})
            case 'freqs'
                freqs = varargin{i+1};
                i = i+2;
            case 'width'
                width = varargin{i+1};
                i = i+2;
            case 'filtfreq'
                filtfreq = varargin{i+1};
                i = i+2;
            case 'filttype'
                filttype = varargin{i+1};
                i = i+2;
            case 'filtorder'
                filtorder = varargin{i+1};
                i = i+2;
            case 'resampledrate'
                resampledrate = varargin{i+1};
                i = i+2;
            case 'downsample'
                dsample = varargin{i+1};
                i = i+2;
            case 'kthresh'
                kthresh = varargin{i+1};
                i = i+2;
            case 'powonly'
                powonly = 1;
                i = i+1;
            case 'usesingles'
                usesingles = 1;
                i = i+1;
            case 'keepk'
                keepk = 1;
                i = i+1;
            otherwise
                error(['Error processing vararg: ' num2str(i)]);
        end
    end
end

% get some parameters
if isempty(freqs)
    freqs = eeganalparams('freqs');
    if isempty(freqs)
        error('EEGTOOLBOX:GETPHASEPOW:NOFREQS','''eeganalparams'' cannot find frequencies');
    end
end
if isempty(width)
    width = eeganalparams('width');
    if isempty(width)
        error('EEGTOOLBOX:GETPHASEPOW:NOFREQS','''eeganalparams'' cannot find width');
    end
end
if isempty(resampledrate)
    % resample all events to the sampling rate of the first event
    resampledrate = GetRateAndFormat(events(1));
end
resampledrate = round(resampledrate);

% convert the durations to samples
buffer = fix((BufferMS)*resampledrate/1000);
% load the eeg
eeg = gete_ms(chan,events,DurationMS+(2*BufferMS),OffsetMS-BufferMS,0,filtfreq,filttype,filtorder,resampledrate);
% keyboard
% see if throw out events with weird kurtosis
if ~isempty(kthresh)
    startsize = size(eeg,1);
    k = kurtosis(eeg(:,buffer+1:end-buffer)');
    goodInd = k<=kthresh;
    %kInd = setdiff(1:size(eeg,1),goodInd);
    kInd = goodInd;
    if keepk
        kInd = k>kthresh;
    else
        eeg = eeg(goodInd,:);
    end
    sizediff = startsize - size(eeg,1);
    if sizediff > 0
        fprintf('Threw out %d events due to kurtosis...\n',sizediff);
    end
else
    if keepk
        kInd = false(size(eeg,1),1);
    else
        kInd = true(size(eeg,1),1);
    end
end

% get the phase and power
[phase,pow] = multiphasevec3(freqs,eeg,resampledrate,width);

% remove the buffer
pow = pow(:,:,buffer+1:end-buffer);
if ~powonly
    phase = phase(:,:,buffer+1:end-buffer);
end


% see if decimate power
if ~isempty(dsample)
    % set the downsampled duration
    dmate = round(resampledrate/dsample);
    dsDur = ceil(size(pow,3)/dmate);
    
    if usesingles
        precision = 'single';
    else
        precision = 'double';
    end
    dpow = zeros(size(pow,1),size(pow,2),dsDur,precision);
    if ~powonly
        dphase = zeros(size(phase,1),size(phase,2),dsDur,precision);
    end
    
    % Must log transform power before decimating
    pow(pow<=0) = eps;
    pow = log10(pow);
    
    % loop and decimate
    fprintf('\nDecimating: %d\n ',size(pow,1));
    for e = 1:size(pow,1)
        fprintf('%d ',e);
        for f = 1:size(pow,2)
            dpow(e,f,:) = decimate(double(pow(e,f,:)),dmate);
            if ~powonly
                % decimate the unwraped phase, then wrap it back
                dphase(e,f,:) = mod(decimate(double(unwrap(phase(e,f,:))),dmate)+pi,2*pi)-pi;
            end
        end
    end
    fprintf('\n');
    
    % replace old pow with new
    pow = dpow;
    clear dpow;
    if ~powonly
        phase = dphase;
        clear dphase;
    end
    
    % convert back to no-log
    pow = 10.^pow;
    
end

if ~powonly
    varargout(1) = {phase};
    varargout(2) = {pow};
    varargout(3) = {kInd};
else
    varargout(1) = {pow};
    varargout(2) = {kInd};
end
