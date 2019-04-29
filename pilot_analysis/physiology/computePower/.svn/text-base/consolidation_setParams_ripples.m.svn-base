function pars = consolidation_setParams_ripples()
% get windowed power for semantic span
pars.froot = '/Volumes/56PUB/UTAH_RAFI/consolidationProject/data/';
pars.eegRootDir = '/Volumes/56PUB/UTAH_RAFI/consolidationProject/data/';
pars.msbuffer = 0;
pars.msmirror = 500;
pars.resamprate = 1000;
pars.bp_flag = 1; %do bipolar referencing
pars.width = 5;
pars.ms_before = [1 1]*1500;
pars.ms_after 	= [1 1]*4500;
pars.ms_before_ret 	= [1 1]*1500;
pars.ms_after_ret 	= [1 1]*4500;
pars.freqs=logspace(log10(70),log10(200),8);

pars.total_points = pars.ms_before + pars.ms_after;
pars.downsample_factor = NaN;%if resamprate is 1000, 20 yields a 50Hz downsampling rate
pars.downsample_points = pars.total_points/pars.downsample_factor;

pars.points_per_win = 100;
pars.points_per_win_ret = 50;
pars.points_per_slide = 50;   %prevent bug propagation


pars.dirWriteOut = '/Volumes/56PUB/UTAH_RAFI/consolidationProject/intermediates/ripples/';
pars.f_name = 'ecog/';

pars.getERP = true;
pars.ERP_dsamsample = 1;
pars.LFP_special = false; %flag for entropy

end