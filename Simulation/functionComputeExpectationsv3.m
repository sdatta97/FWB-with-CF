function [bk, Ck] = ...
    functionComputeExpectationsv3(Hhat_mmW, H_mmW, Hhat_sub6,H_sub6,D,C,nbrOfRealizations,N,N_UE_mmW,N_UE_sub6,K,K_mmW,L,L_mmW,p,gainOverNoise)
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

D_mmW_mmW = zeros(K_mmW,K_mmW,N_UE_mmW,N_UE_mmW);
D_mmW_sub6 = zeros(K_mmW,K-K_mmW,N_UE_mmW,N_UE_sub6);
D_sub6_mmW = zeros(K-K_mmW,K_mmW,N_UE_sub6,N_UE_mmW);
D_sub6_sub6 = zeros(K-K_mmW,K-K_mmW,N_UE_sub6,N_UE_sub6);
Psi_mmW = zeros(K_mmW,N_UE_mmW,N_UE_mmW);
Psi_sub6 = zeros(K-K_mmW,N_UE_sub6,N_UE_sub6);
Psi_mmW_2 = zeros(K_mmW,N_UE_mmW,N_UE_mmW);
Psi_sub6_2 = zeros(K-K_mmW,N_UE_sub6,N_UE_sub6);
%% Compute scaling factors for combining/precoding

%Go through all channel realizations
for n=1:nbrOfRealizations    
    %Go through all APs
    for k = 1:K_mmW
        Psi_mmW (k,:,:) = eye(N_UE_mmW);
        for i = 1:K_mmW
            for l = 1:L
                H = reshape(H_mmW((l-1)*N+1:l*N,n,:,k), [N,N_UE_mmW]);
                Hhat = reshape(Hhat_mmW((l-1)*N+1:l*N,n,:,k), [N,N_UE_mmW]);
                H_int = reshape(H_mmW((l-1)*N+1:l*N,n,:,i), [N,N_UE_mmW]);
                Hhat_int = reshape(Hhat_mmW((l-1)*N+1:l*N,n,:,i), [N,N_UE_mmW]);
                D_mmW_mmW(k,i,:,:) = reshape(D_mmW_mmW(k,i,:,:),[N_UE_mmW,N_UE_mmW]) + (p(k)/(N*L*sum(gainOverNoise(l,:))))*(H'*H_int)./nbrOfRealizations;
            end
%             if (i~=k)
%                 Psi_mmW (k,:,:) = reshape(Psi_mmW(k,:,:),[N_UE_mmW,N_UE_mmW]) + p(k)*reshape(D_mmW_mmW(k,i,:,:),[N_UE_mmW,N_UE_mmW])*reshape(D_mmW_mmW(k,i,:,:),[N_UE_mmW,N_UE_mmW])';
%             end
        end
        for i = 1:K-K_mmW
            for l = 1:L
                H = reshape(H_mmW((l-1)*N+1:l*N,n,:,k), [N,N_UE_mmW]);
                Hhat = reshape(Hhat_mmW((l-1)*N+1:l*N,n,:,k), [N,N_UE_mmW]);
                H_int = reshape(H_sub6((l-1)*N+1:l*N,n,:,i), [N,N_UE_sub6]);
                Hhat_int = reshape(Hhat_sub6((l-1)*N+1:l*N,n,:,i), [N,N_UE_sub6]);
                D_mmW_sub6(k,i,:,:) = reshape(D_mmW_sub6(k,i,:,:),[N_UE_mmW,N_UE_sub6]) + (p(k)/(N*L*sum(gainOverNoise(l,:))))*(H'*H_int)./nbrOfRealizations;
            end
%             Psi_mmW (k,:,:) = reshape(Psi_mmW(k,:,:),[N_UE_mmW,N_UE_mmW]) + p(k)*reshape(D_mmW_sub6(k,i,:,:),[N_UE_mmW,N_UE_sub6])*reshape(D_mmW_sub6(k,i,:,:),[N_UE_mmW,N_UE_sub6])';
        end
%         Psi_mmW_2(k,:,:) = inv(reshape(Psi_mmW(k,:,:),[N_UE_mmW,N_UE_mmW])) + reshape(D_mmW_mmW(k,k,:,:),[N_UE_mmW,N_UE_mmW])*reshape(D_mmW_mmW(k,k,:,:),[N_UE_mmW,N_UE_mmW])';
    end

    for k = 1:K-K_mmW   
        Psi_sub6 (k,:,:) = eye(N_UE_sub6);
        for i = 1:K_mmW
            for l = 1:L
                H = reshape(H_sub6((l-1)*N+1:l*N,n,:,k), [N,N_UE_sub6]);
                Hhat = reshape(Hhat_sub6((l-1)*N+1:l*N,n,:,k), [N,N_UE_sub6]);
                H_int = reshape(H_mmW((l-1)*N+1:l*N,n,:,i), [N,N_UE_mmW]);
                Hhat_int = reshape(Hhat_mmW((l-1)*N+1:l*N,n,:,i), [N,N_UE_mmW]);
                D_sub6_mmW(k,i,:,:) = reshape(D_sub6_mmW(k,i,:,:),[N_UE_sub6,N_UE_mmW]) + (p(k)/(N*L*sum(gainOverNoise(l,:))))*(H'*H_int)./nbrOfRealizations;
            end
%             Psi_sub6 (k,:,:) = reshape(Psi_sub6(k,:,:),[N_UE_sub6,N_UE_sub6]) + p(k)*reshape(D_sub6_mmW(k,i,:,:),[N_UE_sub6,N_UE_mmW])*reshape(D_sub6_mmW(k,i,:,:),[N_UE_sub6,N_UE_mmW])';
        end
        for i = 1:K-K_mmW
            for l = 1:L
                H = reshape(H_sub6((l-1)*N+1:l*N,n,:,k), [N,N_UE_sub6]);
                Hhat = reshape(Hhat_sub6((l-1)*N+1:l*N,n,:,k), [N,N_UE_sub6]);
                H_int = reshape(H_sub6((l-1)*N+1:l*N,n,:,i), [N,N_UE_sub6]);
                Hhat_int = reshape(Hhat_sub6((l-1)*N+1:l*N,n,:,i), [N,N_UE_sub6]);
                D_sub6_sub6(k,i,:,:) = reshape(D_sub6_sub6(k,i,:,:),[N_UE_sub6,N_UE_sub6]) + (p(k)/(N*L*sum(gainOverNoise(l,:))))*(H'*H_int)./nbrOfRealizations;
            end
%             if (i~=k)
%                 Psi_sub6 (k,:,:) = reshape(Psi_sub6(k,:,:),[N_UE_sub6,N_UE_sub6]) + p(k)*reshape(D_sub6_sub6(k,i,:,:),[N_UE_sub6,N_UE_sub6])*reshape(D_sub6_sub6(k,i,:,:),[N_UE_sub6,N_UE_sub6])';
%             end
        end
%         Psi_sub6_2(k,:,:) = inv(reshape(Psi_sub6(k,:,:),[N_UE_sub6,N_UE_sub6])) + reshape(D_sub6_sub6(k,k,:,:),[N_UE_sub6,N_UE_sub6])*reshape(D_sub6_sub6(k,k,:,:),[N_UE_sub6,N_UE_sub6])';
    end            
end
bk = reshape(mean(chgain_arr,1),[L,K]);
Ck = reshape(mean(intgain_arr,1),[L,L,K,K]);
end
