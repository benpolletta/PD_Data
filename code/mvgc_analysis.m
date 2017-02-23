function [varargout] = mvgc_analysis(X,momax,filename,spec_flag)

if nargin < 4 || isempty(spec_flag)
    
    spec_flag = 1;
    
end

[r,c] = size(X);
if r>c
    X = X';
end

[nchans, nobs] = size(X);

if isempty(momax), momax = min(300, nobs - 1); end

% Find optimal model order.
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,'LWR',false);

% figure(1); clf;
% plot_tsdata([AIC BIC]',{'AIC','BIC'},1/1000);
% title('Model order estimation');

morder = moAIC;

% Fit vector autoregressive model.
[A,SIG] = tsdata_to_var(X,morder,'OLS');

assert(~isbad(A),'VAR estimation failed');

% Get autocovariance from VAR.
[G,info] = var_to_autocov(A,SIG,nobs);

switch spec_flag
    
    case 0 % Return Granger.
        
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
        
        varargout{4} = moAIC;
        
        varargout{5} = info;
        
        if ~isempty(filename)
            
            save([filename,'_GC.mat'],'AIC','BIC','moAIC','moBIC','A','SIG','G','info','F','pval','sig')
            
        end
        
    case 1 % Return spectral Granger.
        
        % Get spectral GC from autocovariance.
        if ~isempty(G)
            
            f = autocov_to_spwcgc(G,nobs);
            
            if ~isempty(f)
                
                varargout{1} = permute(f, [3 1 2]);
                
            else
                
                varargout{1} = nan(nobs + 1, nchans, nchans);
                
            end
            
        else
            
            varargout{1} = nan(nobs + 1, nchans, nchans);
            
        end
            
            
        
        varargout{2} = moAIC;
        
        varargout{3} = info;
        
        if ~isempty(filename)
            
            save([filename,'_GC.mat'],'AIC','BIC','moAIC','moBIC','A','SIG','G','info','f')
            
        end
        
    case 2 % Return both Granger and spectral Granger.
        
        % Get Granger from autocovariance.
        F = autocov_to_pwcgc(G);
        
        % varargout{1} = F;
        
        if ~isempty(F)
            
            % Significance test using theoretical null distribution, adjusting for multiple
            % hypotheses.
            pval = mvgc_pval(F,morder,nobs,1,1,1,size(X,1)-2,''); % take careful note of arguments!
            sig  = significance(pval,0.01,'FDR');
            
            % Get spectral GC from autocovariance.
            f = autocov_to_spwcgc(G,nobs);
            
        else
            
            F = nan(size(X,1)); pval = nan(size(X,1)); sig = nan(size(X,1));
            
            f = nan(size(X,1), size(X,1), size(X,2)+1);
            
        end
        
        varargout{1} = F; varargout{2} = pval; varargout{3} = sig;
        
        varargout{4} = f;
        
        varargout{5} = moAIC;
        
        varargout{6} = info;
        
        if ~isempty(filename)
            
            save([filename,'_GC.mat'],'AIC','BIC','moAIC','moBIC','A','SIG','G','info','F','pval','sig','f')
            
        end
        
    case 4 % Return cross-spectrum.
        
        % Get cross-spectrum from autocovariance.
        if ~isempty(G)
            
            [S, ~] = autocov_to_cpsd(G, nobs);
            
            if ~isempty(S)
                
                varargout{1} = permute(S, [3 1 2]);
                
            else
                
                varargout{1} = nan(nobs + 1, nchans, nchans);
                
            end
            
        else
            
            varargout{1} = nan(nobs + 1, nchans, nchans);
            
        end
        
        varargout{2} = moAIC;
        
        varargout{3} = info;
        
        if ~isempty(filename)
            
            save([filename,'_GC.mat'],'AIC','BIC','moAIC','moBIC','A','SIG','G','info','S')
            
        end
        
end

end
