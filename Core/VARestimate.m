function [A,SigmaW,criteria,modal_info] = VARestimate(y,na)

%-- Checking data and splitting into training and validation segments ---
[n,N] = size(y);
N_tra = round(N*0.8);
N_val = N - N_tra;

%-- Building the regression matrix --------------------------------------
Phi = zeros(n*na,N-na);
tau = na+1:N;
for i=1:na
    ind = (1:n) + n*(i-1);
    Phi(ind,:) = y(:,tau-i);
end
Y = y(:,tau);

%-- Estimating AR model -------------------------------------------------
ind_tra = 1:N_tra-na;
ind_val = N_tra-na+(1:N_val);
sss(1) = sum(sum(Y(:,ind_tra).^2));
sss(2) = sum(sum(Y(:,ind_val).^2));

% Computing the AR model parameters and error
A = Y(:,ind_tra)/Phi(:,ind_tra);
Yhat = A*Phi;
err = Y - Yhat;
SigmaW = cov(err');

% Calculating performance criteria
criteria.rss_sss(1) = trace(err(:,ind_tra)*err(:,ind_tra)')/sss(1);
criteria.rss_sss(2) = trace(err(:,ind_val)*err(:,ind_val)')/sss(2);
criteria.lnL = (-1/2)*( N_tra*log(2*pi*det(SigmaW)) + trace( ( SigmaW\err(:,ind_tra) ) * err(:,ind_tra)' ) );
criteria.bic = -2*criteria.lnL + log(N_tra)*numel(A);
criteria.CN = cond(Phi(:,ind_tra));
criteria.spp = N_tra/na;

G = zeros(n*n*na);
for i=na+1:N_tra
    H = kron( Phi(:,i)', eye(n) );
    G = G + (SigmaW\H)'*H;
end
criteria.SigmaTh = (N_tra-na)*pinv(G);

% Calculating modal properties
if nargout > 3
    [V,D] = eig([A; eye((na-1)*n,na*n)]);
    modal_info.phi = V(1:n,:);
    rho = diag(D);
    modal_info.omegaN = abs(log(rho));
    modal_info.zeta = -cos(angle(log(rho)));
end