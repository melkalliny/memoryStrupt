function [elecs, crap] = clean_electrodes(info,pow,bin)

if ~exist('bin','var')
    bin =0;
end

info2 = strrep(info,'PAL1','PAL2');
pow2 = strrep(pow,'PAL1','PAL2');

try
load(info2);
load(pow2);
catch
    warning('unable to find PAL2 session, returning empty matrix')
    elecs = [];
    crap = 1;
    return
end

[~, enc_stim] = filterStruct(evs,'isStim == 1 & strcmp(type,''STUDY_PAIR'') ');
[~, ret_stim] = filterStruct(evs,'isStim == 1 & strcmp(type,''TEST_PROBE'') ');

[~, enc_no_stim] = filterStruct(evs,'isStim == 0 & strcmp(type,''STUDY_PAIR'') ');
[~, ret_no_stim] = filterStruct(evs,'isStim == 0 & strcmp(type,''TEST_PROBE'') ');

if any(~bin)
    error('Shouldn''t be here!')
    d_enc_all = powDataAll(enc_no_stim,:,:);
    s_enc_all = powDataAll(enc_stim,:,:);
    
    d_ret_all = powDataAll(ret_no_stim,:,:);
    s_ret_all = powDataAll(ret_stim,:,:);
    
else
    
    d_enc_all = mean(powDataAll(enc_no_stim,:,:,bin),4);
    s_enc_all = mean(powDataAll(enc_stim,:,:,bin),4);
    
    d_ret_all = mean(powDataAll(ret_no_stim,:,:,bin),4);
    s_ret_all = mean(powDataAll(ret_stim,:,:,bin),4);
end

%%
% el = 95;
% tmp = squeeze(mean(d_enc_all(:,el,:,:),1));
% imagesc(tmp)
% title(num2str(el));
%%


H = NaN(size(powDataAll,2),size(powDataAll,3));
for elec = 1:size(powDataAll,2)
    for f = 1:size(powDataAll,3)
        e_noStim = squeeze(d_enc_all(:,elec,f));
        e_stim   = squeeze(s_enc_all(:,elec,f));
        H(elec,f) = ttest2(e_noStim,e_stim);
    end
end
try
    elecs = ~H;
    crap = 0;
catch eee
    warning('Actually Messed up')
    elecs = ones(size(powDataAll,2),0,0);
   crap = 1;
end