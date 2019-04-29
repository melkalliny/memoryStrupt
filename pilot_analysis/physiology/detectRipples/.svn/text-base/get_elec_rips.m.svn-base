function get_elec_rips(d,p)
% root dir
if ~exist(d.out_dir,'dir'); mkdir(d.out_dir); end

% load events 
load(fullfile(d.task_dir,'events.mat')); 

% load eeg 
load(fullfile(d.out_dir,'eegc.mat'))
load(fullfile(d.out_dir,'mpow.mat')) % just to get windows for ripples

% save hilbert envelope
feeg = buttfilt(eeg,p.rfreq,p.srate,p.type,p.order);
heeg = hilbert(feeg')';
henv = abs(heeg);
zenv = zscore_sess(henv',[ev.session],0);
save(fullfile(d.out_dir,'henv.mat'),'henv','zenv','-v7.3')

% ripple boolean
zenv(zenv<p.thresh)  = 0;
zenv(zenv>=p.thresh) = 1;
zenv = zenv';

% initialize rips
rips  = zeros(size(zenv),'single');
rips_sta  = zeros(size(zenv),'single');
%t_ind = -1*p.off_ms:(-1*p.off_ms+2000);
for i = 1:size(zenv,1)
  series = zenv(i,:);
  rip(i).mid = NaN;
  rip(i).nr  = NaN;
  rip(i).dur = NaN;
  rip(i).isi = NaN;
  if sum(series)<2; continue; end

  % ripple stats
  [sta{i},en{i}] = consecutive_ones(series');
  dur = en{i}-sta{i};
  mid = round(sta{i}+dur/2);
  
  % number duration and isi of ripples
  rip(i).mid = mid;
  rip(i).nr  = sum(dur>p.min_dur);
  rip(i).dur = dur(dur>p.min_dur);
  rip(i).isi = diff(mid(dur>p.min_dur));
  
  % ripple mid
  rips(i,mid) = 1;
  rips_sta(i,sta{i})= 1;
end
m.rips_sta = rips_sta;
m.rips     = rips;
m.rip      = rip;
save(fullfile(d.out_dir,sprintf('rips_%d-%d.mat',p.rfreq)),'-struct','m')



















