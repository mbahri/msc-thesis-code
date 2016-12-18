load('results_yale_welsh_st_sp_0.600000.mat')
mask = ~cellfun(@isempty, results_welsh_st(:));

temp = results_welsh_st(mask);
temp = reshape(temp, 20, 20)';

results_welsh_st = temp
save('results_yale_welsh_st_sp_0.600000.mat', 'results_welsh_st', '-v7.3');