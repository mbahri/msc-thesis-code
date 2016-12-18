% experiment_ten_synth.m
% Executes the experiment on synthetic data with tensor methods.
%
% Georgios Papamakarios
% Imperial College London
% Aug 2014

clear;
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% directory where to write
rootdir = '/vol/bitbucket/gp1813/experiments/synthetic/tensors';
mkdir(rootdir);

% data parameters
dims = {[50 50 50], [100 100 100]};
rank_perc = [0.05 0.1];
data_magn = 1;
err_perc = [0.05 0.1];
err_magn = [100 200];
type = {'tucker_tensor', 'cp_tensor'};

% algorithm parameters
lambda = logspace(-1, -4, 7);
p = [0.1 0.5 1];

for d = 1:length(dims)
    for rp = 1:length(rank_perc)
        for dm = 1:length(data_magn)
            for ep = 1:length(err_perc)
                for em = 1:length(err_magn)
                    for t = 1:length(type)

                        % options
                        opt.dims = dims{d};
                        opt.rank_perc = rank_perc(rp);
                        opt.data_magn = data_magn(dm);
                        opt.err_perc = err_perc(ep);
                        opt.err_magn = err_magn(em);
                        opt.type = type{t};
                        opt.return_error = false;

                        % directory and file to save
                        dir = fullfile(rootdir, opt.type);
                        mkdir(dir)
                        filename = sprintf('%dx%dx%d_rank_%g_data_%g_err_%g_%g', opt.dims, opt.rank_perc, opt.data_magn, opt.err_perc, opt.err_magn);
                        file = fullfile(dir, [filename '.mat']);
                                    
                        % ranks
                        nrank = floor(opt.dims * opt.rank_perc);
                        rank = min(nrank);
                        nrank = {nrank};
                        
                        % create data
                        data = create_synthetic_data(opt);

                        % run experiment
                        fprintf('\n*** %s \n', filename);
                        benchmark_all_ten(data, lambda, p, nrank, rank, false, file);

                    end
                end
            end
        end
    end
end

% save parameters
paramsfile = fullfile(rootdir, 'params.mat');
save(paramsfile, 'dims', 'rank_perc', 'data_magn', 'err_perc', 'err_magn', 'type', 'lambda', 'p');
