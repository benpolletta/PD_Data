function [moAIC,f] = mvgc_analysis(X,momax,filename)

[r,c] = size(X);
if r>c
    X = X';
end

nobs = size(X,2);

% Find optimal model order.
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,'LWR');

morder = moAIC;

% Fit vector autoregressive model.
[A,SIG] = tsdata_to_var(X,morder,'OLS');
    
assert(~isbad(A),'VAR estimation failed');

% Get autocovariance from VAR.
[G,info] = var_to_autocov(A,SIG,nobs);

% Get spectral GC from autocovariance.
f = autocov_to_spwcgc(G,nobs);

if ~isempty(filename)

    save([filename,'_GC.mat'],'AIC','BIC','moAIC','moBIC','A','SIG','G','info','f')
    
end