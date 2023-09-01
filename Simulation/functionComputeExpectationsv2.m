function [signal,signal2, scaling] = ...
    functionComputeExpectationsv2(Hhat_mmW, H_mmW, Hhat_sub6,H_sub6,D,C,nbrOfRealizations,N,N_UE_mmW, N_UE_sub6,K,K_mmW,L,p)
%Compute expectatations that appear in the uplink and downlink SE
%expressions.
%
%INPUT:
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                    with the true channel realizations. The matrix is
%                    organized in the same way as Hhat.
%D                 = DCC matrix with dimension L x K 
%                    where (l,k) is one if AP l serves UE k and zero otherwise
%C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel
%                    estimation error between AP l and UE k,
%                    normalized by noise variance
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Number of UEs per cell
%L                 = Number of APs
%p                 = Vector of UE transmit powers
%
%OUTPUT:
%signal_LP_MMSE    = Matrix with dimension K x K x L
%                    where (i,k,l) is the Monte-Carlo estimation of
%                    expected value of h_{il}^HD_{kl}w_{kl} where w_{kl} is 
%                    LP-MMSE combiner/precoder
%signal2_LP_MMSE   = Matrix with dimension K x K x L
%                    where (i,k,l) is the Monte-Carlo estimation of
%                    expected value of |h_{il}^HD_{kl}w_{kl}|^2 where w_{kl} is 
%                    LP-MMSE combiner/precoder
%scaling_LP_MMSE   = Matrix with dimension L x K
%                    where (l,k) is the Monte-Carlo estimation of
%                    expected value of the norm square of D_{kl}w_{kl}
%                    for LP-MMSE combiner/precoder
%
%This Matlab function was developed to generate simulation results to:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

%Store the N x N identity matrix
eyeN = eye(N);

%Obtain the diagonal matrix with UE transmit powers as its diagonal entries
PowMat = diag(p);

%Scale C by power coefficients
Cp = zeros(size(C));
for k=1:K
    Cp(:,:,:,k) = p(k)*C(:,:,:,k);
end

%Prepare to store simulation results

signal = zeros(K,K,L);
signal2 = zeros(K,K,L);
scaling = zeros(L,K);

chgain_arr = zeros(nbrOfRealizations,L,K);

%% Compute scaling factors for combining/precoding

%Go through all channel realizations
for n=1:nbrOfRealizations    
    %Go through all APs
    for l = 1:L
        for k = 1:K_mmW
            H = reshape(H_mmW((l-1)*N+1:l*N,n_idx,:,k), [N,N_UE_mmW]);
            Hhat = reshape(Hhat_mmW((l-1)*N+1:l*N,n_idx,:,k), [N,N_UE_mmW]);
            [U1, S1, V1] = svd(Hhat');
            [U2, S2, V2] = svd(H');
            chgain_arr(n_idx,l,k) = abs(U1(:,1)'*U2*S2*(V2'*V1(:,1)))^2;
        end
        for k = 1:K-K_mmW
            H = reshape(H_sub6((l-1)*N+1:l*N,n_idx,:,k), [N,N_UE_sub6]);
            Hhat = reshape(Hhat_sub6((l-1)*N+1:l*N,n_idx,:,k), [N,N_UE_sub6]);
            [U1, S1, V1] = svd(Hhat');
            [U2, S2, V2] = svd(H');
            chgain_arr(n_idx,l,k+K_mmW) = abs(U1(:,1)'*U2*S2*(V2'*V1(:,1)))^2;
        end
    end            
end
chgain = reshape(mean(chgain_arr,1),[L,K]);
end
