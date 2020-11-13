function Covariances = LagOneCovSmoother( State, Covariances, System, K, J )

[n,N] = size( State.xtN );

PttmN = zeros(n,n,N);
F = System.F( State.xtt(:,N-1) );
PttmN(:,:,N) = ( eye(n) - K(:,:,N)*System.H )*F*Covariances.Ptt(:,:,N-1);

for i=N:-1:3
    F = System.F( State.xtt(:,i-1) );
    PttmN(:,:,i-1) = Covariances.Ptt(:,:,i-1)*J(:,:,i-2)' + ... 
        J(:,:,i-1)*( PttmN(:,:,i) - F*Covariances.Ptt(:,:,i-1) )*J(:,:,i-2)';
end

Covariances.PttmN = PttmN;