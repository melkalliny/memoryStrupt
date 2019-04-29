
pars = set_params();
global badChan

% [47 49 50 51 52 53 54 57]

for subj = [50 51 52 53 54 57];
    
   subj_dir = sprintf('/Volumes/MK_FRNU/data/eeg/NIH0%d/behavioral/semanticSpan/',subj);
   session_files = dir(subj_dir); sessions = find(contains({session_files.name},'session'));
   montage_dir = sprintf('/Volumes/MK_FRNU/data/eeg/NIH0%d/docs/element_info.csv',subj);
   warning('off')
   montage = readtable(montage_dir);
   montage_dir_bipolar = sprintf('/Volumes/MK_FRNU/data/eeg/NIH0%d/docs/leads_bp.txt',subj);
   fileID = fopen(montage_dir_bipolar,'r');
   bipolar = textscan(fileID,'%s'); bipolar = bipolar{1};
   fclose(fileID);
   pars.electrodes = bipolar;
   
   sessSave = 1; badChan = 0;
   for sessNum = 1:size(sessions,2);
       files = dir(fullfile(subj_dir,session_files(sessions(sessNum)).name));
       istrain = find(contains({files.name},'events.mat'));
       if size(istrain,2) ~= 0;
           events = load(fullfile(subj_dir,session_files(sessions(sessNum)).name,'/events.mat'));
           events = events.events;  
           num_events = size(events,2);
           parfor e = 1:num_events
               try
                   ev = events(e);
                   ev.eegfile = strrep(ev.eegfile,'Shares/FRNU/','MK_FRNU/');
                   compute_power_event_fun_noMir(ev,e,pars,sprintf('NIH0%d',subj),sessSave);
               catch err
                   err.stack(1)
                   err.message
                   fprintf('%d not working\n',e)
                   continue
               end
           end

       end
       sessSave = sessSave+1;
   end

   
   % pull down the bad chans that Julio identifies
   
   % save as subj_session.mat
   fprintf('finished subj %d\n',subj)
   
end

