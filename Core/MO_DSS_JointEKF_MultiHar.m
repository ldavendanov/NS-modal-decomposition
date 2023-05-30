function [Modal,logMarginal] = MO_DSS_JointEKF_MultiHar(y,Orders,HyperPar,Initial)
%--------------------------------------------------------------------------
% Joint EKF estimator for Multiple-Output Diagonal State Space
% representation. This function estimates the 'M' modal components of a
% multicomponent non-stationary signal 'y' using the Joint EKF method.
%
% Created by : David Avendano - September 2020
%    Updated : David Avendano - November 2020
%--------------------------------------------------------------------------

%% Pt 1 : Initial set-up

[d,N] = size(y);                                                            % Size of the signal
p = HyperPar.p;                                                             % Number of frequencies to track
M = numel(Orders)*p;                                                        % Total number of modes to track

% Setting up the state space representation
System.ffun = @(z)ffun(z,Orders,p);
System.F = @(z)stm(z,Orders,p);
System.H = [HyperPar.Psi zeros(d,p)];

%% Pt 2 : Kalman filter for modal estimation

[State,Covariances,logMarginal] = KalmanFilter( y, System, HyperPar, Initial );
State = KalmanSmoother( State, Covariances, System );

fprintf('Log-likelihood %2.4f\n',mean(logMarginal))
    
%% Pt 3 : Packing the results of the modal decomposition

xhat = State.xtN;

Modal.ym = xhat(1:2*M,:);
Modal.theta = xhat(2*M+(1:p),:);
Modal.omega = xhat(2*M+(1:p),:);
Modal.Am = zeros(M,N);
for m=1:M
    Modal.Am(m,:) = sqrt( xhat(2*m-1,:).^2 + xhat(2*m,:).^2 );
end

% % Making sure that the modes are normalized
% C = zeros(M,1);
% for m=1:M
%     C(m) = std(Modal.ym(2*m,:));
%     Modal.ym((1:2)+2*(m-1),:) = Modal.ym((1:2)+2*(m-1),:)/C(m);
%     Modal.Am(m,:) = Modal.Am(m,:)/C(m);
%     HyperPar.Psi(:,(1:2)+2*(m-1)) = C(m)*HyperPar.Psi(:,(1:2)+2*(m-1));
% end
% 
%--------------------------------------------------------------------------
function z_new = ffun(z_old,ord,p)

n = length(z_old)-p;
M = dffun_dz(z_old(n+(1:p)),ord,p);
z_new = z_old;
z_new(1:n) = M*z_old(1:n);

%--------------------------------------------------------------------------
function F = stm(z,ord,p)

n = length(z)-p;
M = dffun_dz(z(n+(1:p)),ord,p);
Z = dffun_dtheta(z(1:n),z(n+(1:p)),ord,p);
F = [M Z; zeros(p,n) eye(p)];

%--------------------------------------------------------------------------
function M = dffun_dz(theta,ord,p)

n = length(ord);
M = zeros(2*n*p);

for j=1:p
    for k=1:n
        ind = (1:2)  +2*(k-1) + 2*n*(j-1);
        M(ind,ind) = [ cos(ord(k)*theta(j)) sin(ord(k)*theta(j))
                      -sin(ord(k)*theta(j)) cos(ord(k)*theta(j))];
    end
end

%--------------------------------------------------------------------------
function Z = dffun_dtheta(z,theta,ord,p)

n = length(ord);
Z = zeros(2*n*p,p);

for j=1:p
    for k=1:n
        ind = (1:2) + 2*(k-1) + 2*n*(j-1);
        Z(ind,j) = ord(k)*[-sin(ord(k)*theta(j))  cos(ord(k)*theta(j));
                           -cos(ord(k)*theta(j)) -sin(ord(k)*theta(j))]*z(ind);
    end
end

% 
% 
% %--------------------------------------------------------------------------
% function z_new = ffun(z_old,ord,omega0)
% 
% n = length(z_old)-1;
% M = dffun_dz(z_old(end),ord,omega0);
% z_new = z_old;
% z_new(1:n) = M*z_old(1:n);
% 
% %--------------------------------------------------------------------------
% function F = stm(z,ord,omega0)
% 
% n = length(z)-1;
% M = dffun_dz(z(end),ord,omega0);
% Z = dffun_dtheta(z(1:n),z(end),ord,omega0);
% F = [M Z; zeros(1,n) 1];
% 
% %--------------------------------------------------------------------------
% function M = dffun_dz(theta,ord,omega0)
% 
% n = length(ord);
% M = zeros(2*n);
% theta = theta + omega0;
% for k=1:n
%     ind = (1:2)+2*(k-1);
%     M(ind,ind) = [ cos(ord(k)*theta) sin(ord(k)*theta)
%                   -sin(ord(k)*theta) cos(ord(k)*theta)];
% end
% 
% %--------------------------------------------------------------------------
% function Z = dffun_dtheta(z,theta,ord,omega0)
% 
% n = length(ord);
% Z = zeros(2*n,1);
% theta = theta + omega0;
% for k=1:n
%     ind = (1:2)+2*(k-1);
%     Z(ind) = ord(k)*[-sin(ord(k)*theta)  cos(ord(k)*theta);
%                      -cos(ord(k)*theta) -sin(ord(k)*theta)]*z(ind);
% end
