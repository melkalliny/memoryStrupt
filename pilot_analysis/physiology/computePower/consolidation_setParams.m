function pars = consolidation_setParams()
% get windowed power for semantic span
pars.froot = '/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/';
pars.eegRootDir = '/Volumes/Shares/FRNU/dataWorking/consolidationProject/data/';
pars.msbuffer = 0;
pars.msmirror = 500;
pars.resamprate = 1000;
pars.bp_flag = 1; %do bipolar referencing
pars.width = 5;
pars.ms_before = [1 1]*1500;
pars.ms_after 	= [1 1]*4500;
pars.ms_before_ret 	= [1 1]*4500;
pars.ms_after_ret 	= [1 1]*1500;
%pars.freqs=logspace(log10(3),log10(200),70);
pars.freqs=cat(2,logspace(log10(4),log10(12),10),logspace(log10(80),log10(150),10));

pars.total_points = pars.ms_before + pars.ms_after;
pars.downsample_factor = NaN;%if resamprate is 1000, 20 yields a 50Hz downsampling rate
pars.downsample_points = pars.total_points/pars.downsample_factor;

pars.points_per_win = 100;
pars.points_per_win_ret = 100;
pars.points_per_slide = 10;   %prevent bug propagation

pars.BipolarOrMonopolar = 'monopolar'; 
pars.dirWriteOut = '/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/power_LowAndHigh_monopolar_allPoints/';
pars.f_name = 'ecog/';

pars.getERP = true;
pars.ERP_dsamsample = 1;
pars.LFP_special = false; %flag for entropy

end