function [State,Covariances,J] = KalmanSmoother( State, Covariances, System )

% Unpacking input
xtt = State.xtt;            xttm = State.xttm;
Ptt = Covariances.Ptt;      Pttm = Covariances.Pttm;
[n,N] = size(xtt);

% Initializing
xtN = xtt;
PtN = Ptt;
J = zeros(n,n,N);

for i=N:-1:2
    F = System.F( xtt(:,i-1) );
    J(:,:,i-1) = Ptt(:,:,i-1)*F'/Pttm(:,:,i);
    xtN(:,i-1) = xtt(:,i-1) + J(:,:,i-1)*( xtN(:,i) - xttm(:,i) );
    PtN(:,:,i-1) = Ptt(:,:,i-1) + J(:,:,i-1)*( PtN(:,:,i) - Pttm(:,:,i) )*J(:,:,i-1)';
end

% Packing results
State.xtN = xtN;
Covariances.PtN = PtN;