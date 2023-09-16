function [bk, Ck] = ...
    functionComputeExpectationsv4(Hhat_mmW, H_mmW, Hhat_sub6,H_sub6,D,C,nbrOfRealizations,N,N_UE_mmW,N_UE_sub6,K,K_mmW,L,L_mmW,p,gainOverNoise)
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
eta_eq_power = 1./(N*L*sum(gainOverNoise,2));
bk = gainOverNoise;
Ck = zeros(L,L,K,K);
for k = 1:K_mmW
    for i = 1:K_mmW
        Ck (:,:,k,i) = (N_UE_mmW/N)*gainOverNoise(:,k)*gainOverNoise(:,i)';
        if (i==k)
            for l = 1:L
                for ll = 1:L
                    if (ll == l)
                        Ck (l,ll,k,i) = Ck (l,ll,k,i) + (gainOverNoise(l,k))^2;
                    else 
                        Ck (l,ll,k,i) = Ck (l,ll,k,i) + 2*gainOverNoise(l,k)*gainOverNoise(ll,k);
                    end

                end
            end
        end
    end
    for i = 1:K-K_mmW
        Ck (:,:,k,i+K_mmW) = (N_UE_sub6/N)*gainOverNoise(:,k)*gainOverNoise(:,i+K_mmW)';
    end
end
for k = 1:K-K_mmW
    for i = 1:K_mmW
        Ck (:,:,k+K_mmW,i) = (N_UE_mmW/N)*gainOverNoise(:,k+K_mmW)*gainOverNoise(:,i)';
    end
    for i = 1:K-K_mmW
        Ck (:,:,k+K_mmW,i+K_mmW) = (N_UE_sub6/N)*gainOverNoise(:,k+K_mmW)*gainOverNoise(:,i+K_mmW)';
        if (i==k)
            for l = 1:L
                for ll = 1:L
                    if (ll == l)
                        Ck (l,ll,k+K_mmW,i+K_mmW) = Ck (l,ll,k+K_mmW,i+K_mmW) + (gainOverNoise(l,k+K_mmW))^2;
                    else 
                        Ck (l,ll,k+K_mmW,i+K_mmW) = Ck (l,ll,k+K_mmW,i+K_mmW) + 2*gainOverNoise(l,k+K_mmW)*gainOverNoise(ll,k+K_mmW);
                    end
                end
            end
        end
    end
end
end
