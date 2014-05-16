function [moAIC,info,varargout] = mvgc_analysis(X,momax,filename,spec_flag)

if nargin < 4 || isempty(spec_flag)
    
    spec_flag = 1;
    
end

[r,c] = size(X);
if r>c
    X = X';
end

nobs = size(X,2);

% Find optimal model order.
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,'LWR',false);

morder = moAIC;

% Fit vector autoregressive model.
[A,SIG] = tsdata_to_var(X,morder,'OLS');
    
assert(~isbad(A),'VAR estimation failed');

% Get autocovariance from VAR.
[G,info] = var_to_autocov(A,SIG,nobs);

if spec_flag == 0
    
    % Get Granger from autocovariance.
    F = autocov_to_pwcgc(G);
    
    varargout{1} = F;
    
    % Significance test using theoretical null distribution, adjusting for multiple
    % hypotheses.
    
    if ~isempty(F)
        
        pval = mvgc_pval(F,morder,nobs,1,1,1,size(X,1)-2,''); % take careful note of arguments!
        sig  = significance(pval,0.01,'FDR');
        
    else
        
        pval = nan; sig = nan;
        
    end
        
    varargout{2} = pval; varargout{3} = sig;
    
    if ~isempty(filename)
        
        save([filename,'_GC.mat'],'AIC','BIC','moAIC','moBIC','A','SIG','G','info','F','pval','sig')
        
    end
    
elseif spec_flag == 1
    
    % Get spectral GC from autocovariance.
    f = autocov_to_spwcgc(G,nobs);
    
    varargout{1} = f;
    
    if ~isempty(filename)
        
        save([filename,'_GC.mat'],'AIC','BIC','moAIC','moBIC','A','SIG','G','info','f')
        
    end
    
end