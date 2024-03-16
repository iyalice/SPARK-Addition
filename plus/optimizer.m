load GAPDH.mat
X = double([X.x, X.y]);
y = double(y);

%% Load GPML
addpath(genpath('plus\gpml'))

%% Set up model.
meanfunc = {@meanZero};
hyp.mean = [];

covfunc1 = {@covSum,{{@covMask,{[1 0],@covSEiso}},{@covMask,{[0 1],@covSEiso}}}}; 
hyp.cov = [0 0 0 0];

likfunc = @likGauss;
hyp.lik = log(var(y)/10);

%% Repeat...
[hyp_opt, nlls] = minimize(hyp, @gp, -int32(100 * 3 / 3), @infExact, meanfunc, covfunc1, likfunc, X, y);
best_nll = nlls(end);

fid1 = fopen('./plus/GAPDH_1.in','w');
fprintf(fid1,'exp(-(X[i,1]-X[j,1])^2/(2*exp(%f)^2))*exp(%f)^2+exp(-(X[i,2]-X[j,2])^2/(2*exp(%f)^2))*exp(%f)^2',hyp_opt.cov);

covfunc2 = {@covProd,{{@covMask,{[1 0],@covSEiso}},{@covMask,{[0 1],@covSEiso}}}}; 
hyp2.mean = [];
hyp2.cov = [0 0 0 0];
hyp2.lik = log(var(y)/10);
[hyp2_opt, nlls] = minimize(hyp2, @gp, -int32(100 * 3 / 3), @infExact, meanfunc, covfunc2, likfunc, X, y);
best_nll2 = nlls(end);


fid2 = fopen('./plus/GAPDH_2.in','w');
fprintf(fid2,'exp(-(X[i,1]-X[j,1])^2/(2*exp(%f)^2))*exp(%f)^2*exp(-(X[i,2]-X[j,2])^2/(2*exp(%f)^2))*exp(%f)^2',hyp2_opt.cov);
