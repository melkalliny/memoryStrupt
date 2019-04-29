
writeBiowulf = 1;
cmd_num = 1;

swarm_fname = 'calculate_swarm.sh';
big_bash_fname = 'calculate_big_bash.sh';

biowulf_swarm_script_dir = '/home/elkallinymm/consolidationProject/biowulf_packaging/calculate/scripts';
biowulf_swarm_fpath = [biowulf_swarm_script_dir '/' swarm_fname];
biowulf_big_bash_fpath = [biowulf_swarm_script_dir '/' big_bash_fname];

local_swarm_script_dir = '/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/code/biowulf_packaging/scripts/';
script_big_bash = [ local_swarm_script_dir big_bash_fname];
script_swarm = [local_swarm_script_dir swarm_fname];

big_bash_fid = fopen(script_big_bash, 'w');
swarm_fid = fopen(script_swarm, 'w');

fprintf(swarm_fid, [ 'mkdir ' biowulf_swarm_script_dir '/log_dump\n' ]);
fprintf(swarm_fid, 'swarm -g 50 -b 1 -t 1 --time 12:00:00 --gres=lscratch:15 --merge-output --logdir %s -f  %s\n', [ biowulf_swarm_script_dir '/log_dump' ],  biowulf_big_bash_fpath);

windowLength = [100 250 500 1000 2000];
windowShift = [10 25 50 100 200];
for windowOptions = 1:3
    subjectNames = {'NIH043', 'NIH050', 'NIH051', 'NIH052', 'NIH053', 'NIH055'};
    params = [1 2 3 4 5 6 1 2 3 4 5 6; 2 2 2 2 2 2 1 1 1 1 1 1];
    
    correlations_results = cell(1,12);
    subjectsToRun = [1:6];
    for doParam = [1:12]
        strings = {'PROBE_START','RESPONSE'};
        params = [1 2 3 4 5 6 1 2 3 4 5 6; 2 2 2 2 2 2 1 1 1 1 1 1];
        stringEvent = strings{params(2,doParam)};
        
        if writeBiowulf
                        
            biowulf_matlab_script_dir = '/home/elkallinymm/consolidationProject/biowulf_packaging/calculate';
            biowulf_matlab_script_bash_name = 'run_consolidation_calculateSimilarity_LowAndHigh_biowulf_swarm.sh';
            
            lib_string = 'tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v94.tar.gz;';
            
            cmd_string = sprintf('cd %s; ./%s /lscratch/$SLURM_JOB_ID/v94 ', biowulf_matlab_script_dir, biowulf_matlab_script_bash_name);
            cmd_string = sprintf('%s %s %s %s %s', cmd_string, subjectNames{params(2,doParam)}, num2str(subjectsToRun(params(1,doParam))), stringEvent, num2str(windowOptions));
            
            
            biowulf_cmd_fpath = [ biowulf_swarm_script_dir '/' num2str(cmd_num) '.sh' ];
            biowulf_log_fpath = [ biowulf_swarm_script_dir '/' num2str(cmd_num) '.log' ];
            
            local_cmd_fpath = [ local_swarm_script_dir '/' num2str(cmd_num) '.sh' ];
            
            local_cmd_fid = fopen(local_cmd_fpath, 'w');
            
            fprintf(local_cmd_fid, '#!/bin/bash\n');
            fprintf(local_cmd_fid, '%s\n', lib_string);
            fprintf(local_cmd_fid, '%s &> %s\n', cmd_string, biowulf_log_fpath);
            
            fclose(local_cmd_fid);
            
            fprintf(big_bash_fid, '%s\n', ['bash ' biowulf_cmd_fpath]);
            
            cmd_num = cmd_num + 1;
            
        else
            %[correlations_results{1,doParam}] = consolidation_calculateSimilarity_LowAndHigh(subjectNames, num2str(subjectsToRun(params(1,doParam))),stringEvent, windowLength(windowOptions), windowShift(windowOptions));
        end
    end
    %save(sprintf('/Volumes/FRNU_SVN/consolidation_SVN/consolidationProject/intermediates/correlations_windowSizes/windowSize%d.mat',windowOptions),'correlations_results')
    
end

fclose(big_bash_fid);
fclose(swarm_fid);