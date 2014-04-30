
%% Parameters

ntrials   = 1;     % number of trials
% nobs      = 1000;   % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

momax     = 300;     % maximum model order for model order estimation

% acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 1000;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

%% Load data.

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

basetimes = [300 1200 1800 600 1800];

infusetimes = [390 240 300 450 510];

load('BetaTimes')

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    load([folder,'/',prefix,'_beta_data.mat']);
    
    X = beta_data';
    
    [nvars, nobs] = size(X);
    
    acmaxlags = size(X,2);
    
    %% Model order estimation
    
    % Calculate information criteria up to specified maximum model order.
    
    if isempty(dir([folder,'/',prefix,'_beta_data_infocrit.mat']))
        
        ptic('\n*** tsdata_to_infocrit\n');
        [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
        ptoc('*** tsdata_to_infocrit took ');
        
        save([folder,'/',prefix,'_beta_data_infocrit.mat'],'AIC','BIC','moAIC','moBIC')
    
    else
        
        load([folder,'/',prefix,'_beta_data_infocrit.mat'])
        
    end
    
    % Plot information criteria.
    
    figure(1); clf;
    plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    title('Model order estimation');
    
    save_as_pdf(gcf,[folder,'/',prefix,'_beta_data_infocrit']);
    
    fprintf('\nbest model order (AIC) = %d\n',moAIC);
    fprintf('best model order (BIC) = %d\n',moBIC);
    
    % Select model order.
    
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
    
    %% VAR model estimation
    
    % Estimate VAR model of selected order from data.
    
    if isempty(dir([folder,'/',prefix,'_beta_data_VAR.mat']))
        
        ptic('\n*** tsdata_to_var... ');
        [A,SIG] = tsdata_to_var(X,morder,regmode);
        % [F,A,SIG,E] = GCCA_tsdata_to_pwcgc(X',morder,regmode);
        ptoc;
        
        save([folder,'/',prefix,'_beta_data_VAR.mat'],'A','SIG')
        
    else
        
        load([folder,'/',prefix,'_beta_data_VAR.mat'])
        
    end
    
    % Check for failed regression
    
    assert(~isbad(A),'VAR estimation failed');
    
    % NOTE: at this point we have a model and are finished with the data! - all
    % subsequent calculations work from the estimated VAR parameters A and SIG.
    
    %% Autocovariance calculation
    
    % The autocovariance sequence drives many Granger causality calculations (see
    % next section). Now we calculate the autocovariance sequence G according to the
    % VAR model, to as many lags as it takes to decay to below the numerical
    % tolerance level, or to acmaxlags lags if specified (i.e. non-empty).
    
    if isempty(dir([folder,'/',prefix,'_beta_data_AC.mat']))
        
        ptic('*** var_to_autocov... ');
        [G,info] = var_to_autocov(A,SIG,acmaxlags);
        ptoc;
        
        save([folder,'/',prefix,'_beta_data_AC.mat'],'G','info')
        
    else
        
        load([folder,'/',prefix,'_beta_data_AC.mat'])
        
    end
    
    % The above routine does a LOT of error checking and issues useful diagnostics.
    % If there are problems with your data (e.g. non-stationarity, colinearity,
    % etc.) there's a good chance it'll show up at this point - and the diagnostics
    % may supply useful information as to what went wrong. It is thus essential to
    % report and check for errors here.
    
    var_info(info,true); % report results (and bail out on error)
    
    %% Granger causality calculation: time domain
    
    if isempty(dir([folder,'/',prefix,'_beta_data_PWGC.mat']))
        
        % Calculate time-domain pairwise-conditional causalities - this just requires
        % the autocovariance sequence.
        
        ptic('*** autocov_to_pwcgc... ');
        F = autocov_to_pwcgc(G);
        ptoc;
        
        % Significance test using theoretical null distribution, adjusting for multiple
        % hypotheses.
        
        pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
        sig  = significance(pval,alpha,mhtc);
        
        save([folder,'/',prefix,'_beta_data_PWGC.mat'],'F','pval','sig')
        
    else
        
        load([folder,'/',prefix,'_beta_data_PWGC.mat'])
        
    end
    
    % Check for failed GC calculation
    
    assert(~isbad(F,false),'GC calculation failed');
    
    % Plot time-domain causal graph, p-values and significance.
    
    figure(2); clf;
    subplot(1,3,1);
    plot_pw(F);
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(pval);
    title('p-values');
    subplot(1,3,3);
    plot_pw(sig);
    title(['Significant at p = ' num2str(alpha)])
    
    save_as_pdf(gcf,[folder,'/',prefix,'_beta_data_PWGC'])
    
    % For good measure we calculate Seth's causal density (cd) measure - the mean
    % pairwise-conditional causality. We don't have a theoretical sampling
    % distribution for this.
    
    cd = mean(F(~isnan(F)));
    
    fprintf('\ncausal density = %f\n',cd);
    
    %% Granger causality calculation: frequency domain
    
    % Calculate spectral pairwise-conditional causalities at given frequency
    % resolution - again, this only requires the autocovariance sequence.
    
    if isempty(dir([folder,'/',prefix,'_beta_data_GCspec.mat']))
        
        ptic('\n*** autocov_to_spwcgc... ');
        f = autocov_to_spwcgc(G,fres);
        ptoc;
        
        save([folder,'/',prefix,'_beta_data_GCspec.mat'])
        
    else
        
        load([folder,'/',prefix,'_beta_data_GCspec.mat'])
        
    end
    
    % Check for failed spectral GC calculation
    
    assert(~isbad(f,false),'spectral GC calculation failed');
    
    % Plot spectral causal graph.
    
    figure(3); clf;
    plot_spw_1axis(f,fs,[10 40]);
    title([folder,' Spectral Granger Causality'])
    
    save_as_pdf(gcf,[folder,'/',prefix,'_beta_data_GCspec'])
    
    %% Granger causality calculation: frequency domain -> time-domain  (<mvgc_schema.html#3 |A15|>)
    
    % Check that spectral causalities average (integrate) to time-domain
    % causalities, as they should according to theory.
    
    fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
    Fint = smvgc_to_mvgc(f); % integrate spectral MVGCs
    mad = maxabs(F-Fint);
    madthreshold = 1e-5;
    if mad < madthreshold
        fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
    else
        fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
    end
    
end
