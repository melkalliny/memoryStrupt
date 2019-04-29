function EEG=gete(channel,events,duration,offset,buffer,filtfreq,filttype,filtorder,resampledrate)
%GETE - Get EEG event data.
% 
% Returns data from an eeg file.  User specifies the channel,
% duration, and offset along with an event.  The event struct MUST
% contain both 'eegfile' and 'eegoffset' members.
%
% If the duration is set to 0, then it will return a cell array of
% the entire data from each unique file in the events structure.
% It will ignore the offset and buffer variables.
%
% FUNCTION:
%   EEG=gete(channel,events,duration,offset,buffer,filtfreq,filttype,filtorder,resampledrate)
%
% INPUT ARGS:
%   channel = 'G28';            % channel name (string) or cell string (for bipolar)
%   events = events(8:12);  % event struct to extract [eegfile eegoffset]
%   duration = 512;         % signal time length in samples (of the
%                           % original duration of the data)
%   offset = 0;             % offset at which to start in samples
%                           % (in units of the original data)
%   buffer = 256;           % buffer (needed for filtering)
%                           %   default is 0
%   filtfreq = [59.5 60.5]; % Filter freq (depends on type, see buttfilt)
%                           %   default is []
%   filttype = 'stop';      % Filter type (see buttfilt)
%   filtorder = 1;          % Filter order (see buttfilt)
%   resampledrate = 256     % (optional) - resampled the data
%   remap_sh                % (optional) - if 1, if given channel does not exist, check for its _sh and return
%
% OUTPUT ARGS:
%   EEG - The data from the file
%

% 2016/10/24 - MST: allow channel names instead of numbers
% 2010/7/2 - JRM: return nans if event offset is past the end of the file
% 2006/8/4 - MvV: fixed bug in case of fractionated durations
% 2006/1/20 - MvV: added in resampling, similar to gete_ms.m
% 2004/3/5 - PBS: Now reads in entire unique file if duration is 0
% 2003/12/9 - PBS: Added which event number to read warning

% check the arg
if nargin < 9
  resampledrate = [];
if nargin < 8
  filtorder = 1;
  if nargin < 7
    filttype = 'stop';
    if nargin < 6
      filtfreq = [];
      if nargin < 5
	buffer = 0;
	if nargin<4 
	  offset=0; 
	end; 
      end
    end
  end
end
end

% dereference char cells
if iscellstr(channel) && isscalar(channel)
    channel = char(channel);
end

% get data info
[samplerate,nBytes,dataformat,gain] = GetRateAndFormat(events(1));
samplerate = round(samplerate);

if isempty(resampledrate)
  resampledrate = samplerate;
end
resampledrate = round(resampledrate);

final_duration = fix(duration*resampledrate/samplerate);
final_buffer = fix(buffer*resampledrate/samplerate);


% see if getting for each event or all unique files
if duration == 0
  % getting for all unique files
  uFiles = unique(getStructField(events,'eegfile',''));
  if isempty(uFiles{1})
      uFiles = uFiles(2:end);
  end
  
  EEG = cell(1,length(uFiles));
  
  % loop over files
  for f = 1:length(uFiles)
        % set the channel filename
        if iscellstr(channel) && numel(channel) == 2
            % Bipolar
            eegfname=fullfile(uFiles{f}, sprintf('%s-%s', channel{:}));
        else
            % Monopolar
           eegfname=fullfile(uFiles{f}, channel);
        end
        
        eegfile=fopen(eegfname,'r','l'); % NOTE: the 'l' means that it came from a PC!
        
        
        % tell if not open
        while eegfile==-1
            % did not open
          
            if strfound(eegfname, '_sh')
                eegfname_unshift = strrep(eegfname, '_sh', '');
                eegfile = fopen(eegfname_unshift,'r','l');
                warning('EEG file not found, ** using original (without _sh) instead! ** : %s', eegfile);
                
            else
                error('ERROR: EEG File not found: %s.\n',eegfname);
            end
          
          
          
        end
  
        % read the entire file
        EEG{f}= fread(eegfile,inf,dataformat)';
        
        % close the file
        try 
            fclose(eegfile);
        catch, end
    
        % see if filter the data
        if ~isempty(filtfreq)
            EEG{f}=buttfilt(EEG{f},filtfreq,samplerate,filttype,filtorder);
        end
  
	% see if resample
	if resampledrate ~= samplerate
	  % do the resample
	  EEG{f} = resample(EEG{f},round(resampledrate),round(samplerate));
	  %readbytes = dsamp(readbytes,round(samplerate),round(resampledRate));
	end	
	
        % apply the gain
        EEG{f} = EEG{f}.*gain;
  end
  
else 
  % getting for each event
  % allocate space
  EEG = zeros(length(events),final_duration+(2*final_buffer));
  
  for e = 1:length(events)


    % set the channel filename
    format = ['%s/' turnary(isnumeric(channel), '%03i', '%s')];
    eegfname = sprintf(format, events(e).eegfile, channel); % '%s/%s' or '%s/%03i'
    eegfile = fopen(eegfname,'r','l'); % NOTE: the 'l' means that it came from a PC!
    
    % tell if not open
    if eegfile==-1
      % did not open
      error('ERROR: EEG File not found for event(%d): %s.\n',e,events(e).eegfile);
    end
    
    % read the eeg data
    thetime=offset+events(e).eegoffset-buffer;
    status = fseek(eegfile,nBytes*thetime,-1);    
    
    if status == 0
        readbytes = fread(eegfile,duration+(2*buffer),dataformat)';
    elseif status == -1
        warning('EEGTOOLBOX:GETE:NODATA','%s: eeg data for event %d were not found',eegfname,e);
        readbytes = nan(1,duration+(2*buffer));
    end

    if length(readbytes)~=fix(duration+2*buffer)
      warning('EEGTOOLBOX:GETE:INCOMPLETEDATA','%s: only %d of %d samples read for event %d -- appending nans',eegfname,length(readbytes),size(EEG,2),e);
      readbytes = [readbytes nan(1,size(EEG,2)-length(readbytes))]; %#ok<AGROW>
    end
        
    % close the file
    fclose(eegfile);
    
    % see if filter the data
    if ~isempty(filtfreq)
      readbytes=buttfilt(readbytes,filtfreq,samplerate,filttype,filtorder);
    end
    % see if resample
    if resampledrate ~= samplerate
      % do the resample
      readbytes = resample(readbytes,round(resampledrate),round(samplerate));
      %readbytes = dsamp(readbytes,round(samplerate),round(resampledrate));
    end
    
    EEG(e,:)=readbytes;  
    
  end
  
  EEG = EEG(:,final_buffer+1:end-final_buffer);
  
  % apply the gain
  EEG = EEG.*gain;
end




