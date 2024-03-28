
%% Load GPML
addpath(genpath('gpml'))
gene_table = readtable('./BC_pvals.csv');
for k = 1:107
    gene = cell2mat(gene_table{k,"gene"});
    path = sprintf('./%s.mat' ,gene);
    
    load(path)
    X = double([X.x,X.y]);
    y = double(y);
    %% Set up model.
    meanfunc = {@meanZero};
    hyp.mean = [];
    likfunc = @likGauss;
    covfunc = @covSEisoU; 
    hyp.mean = [];
    hyp.cov = 0;
    hyp.lik = log(var(y)/10);
    [hyp_opt, nlls] = minimize(hyp, @gp, -int32(100 * 3 / 3), @infExact, meanfunc, covfunc, likfunc, X, y);
    best_nll = nlls(end);
    
    
    fid = fopen("./new_result_2.csv",'a');
    fprintf(fid,'%s,%f\n',gene,hyp_opt.cov);
end