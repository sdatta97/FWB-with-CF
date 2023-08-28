function [signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
    functionComputeExpectations(Hhat_mmW,Hhat_sub6,H_mmW,H_sub6,D,C,nbrOfRealizations,N,N_UE_mmW,N_UE_sub6,K_mmW,K,L_mmW,L,p)
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

signal_LP_MMSE = zeros(K,K,L);
signal2_LP_MMSE = zeros(K,K,L);
scaling_LP_MMSE = zeros(L,K);

%% Compute scaling factors for combining/precoding

%Go through all channel realizations
for n=1:nbrOfRealizations    
    %Go through all APs
    for l = 1:L
        %Extract channel realizations from all UEs to AP l
        % Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);
        Hallj_mmW = reshape(H_mmW(1+(l-1)*N:l*N,n,:,:),[N N_UE_mmW K_mmW]);
        Hallj_sub6 = reshape(H_sub6(1+(l-1)*N:l*N,n,:,:),[N N_UE_sub6 K-K_mmW]);

        %Extract channel estimate realizations from all UEs to AP l
        % Hhatallj = reshape(Hhat(1+(l-1)*N:l*N,n,:),[N K]);
        Hhatallj_mmW = reshape(Hhat_mmW(1+(l-1)*N:l*N,n,:,:),[N N_UE_mmW K_mmW]);        
        Hhatallj_sub6 = reshape(Hhat_sub6(1+(l-1)*N:l*N,n,:,:),[N N_UE_sub6 K-K_mmW]);        

        %Extract which UEs are served by AP l
        % servedUEs = find(D(l,:)==1);
        servedUEs_mmW = find(D(l,1:K_mmW)==1);
        servedUEs_sub6 = find(D(l,1+K_mmW:K)==1);
        %Obtain the statistical matrices used for
        %computing partial combining/precoding schemes
        % Cpserved = reshape(sum(Cp(:,:,l,servedUEs),4),[N N]);
        % Pserved = PowMat(servedUEs,servedUEs);
        Cpserved_mmW = reshape(sum(Cp(:,:,l,servedUEs_mmW),4),[N N]);
        Pserved_mmW = PowMat(servedUEs_mmW,servedUEs_mmW);
        Cpserved_sub6 = reshape(sum(Cp(:,:,l,K_mmW+servedUEs_sub6),4),[N N]);
        Pserved_sub6 = PowMat(K_mmW+servedUEs_sub6,K_mmW+servedUEs_sub6);

        %Compute MR combining scaled by square root of transmit powers
        % Vp_MR = Hhatallj(:,servedUEs)*sqrt(Pserved);
        Vp_MR_mmW  = Hhatallj_mmW(:,:,servedUEs_mmW)*sqrt(Pserved_mmW);
        Vp_MR_sub6 = Hhatallj_sub6(:,:,servedUEs_sub6)*sqrt(Pserved_sub6);
        %Compute LP-MMSE combining
        % V_LP_MMSE = (((Vp_MR*Vp_MR')+Cpserved+eyeN)\Vp_MR)*sqrt(Pserved);
        V_LP_MMSE_mmW = (((Vp_MR_mmW*Vp_MR_mmW')+Cpserved_mmW+eyeN)\Vp_MR_mmW)*sqrt(Pserved_mmW);
        V_LP_MMSE_sub6 = (((Vp_MR_sub6*Vp_MR_sub6')+Cpserved_sub6+eyeN)\Vp_MR_sub6)*sqrt(Pserved_sub6);
        
        %Go through all UEs served by the AP
        % for ind = 1:length(servedUEs)
        % 
        %     %Extract UE index
        %     % k = servedUEs(ind);
        %     k = servedUEs(ind);
        % 
        %     %Normalize LP-MMSE precoding
        %     w = V_LP_MMSE(:,ind);
        % 
        %     %Compute realizations of the terms inside the expectations
        %     %of the signal and interference terms in the SE expressions and
        %     %update Monte-Carlo estimates 
        %     signal2_LP_MMSE(:,k,l) = signal2_LP_MMSE(:,k,l) + abs(Hallj'*w).^2/nbrOfRealizations;            
        %     signal_LP_MMSE(:,k,l) = signal_LP_MMSE(:,k,l) + Hallj'*w/nbrOfRealizations;           
        %     scaling_LP_MMSE(l,k) = scaling_LP_MMSE(l,k) + sum(abs(w).^2,1)/nbrOfRealizations;            
        % end
        for ind = 1:length(servedUEs_mmW)
            
            %Extract UE index
            % k = servedUEs(ind);
            k = servedUEs_mmW(ind);
            
            %Normalize LP-MMSE precoding
            w = V_LP_MMSE_mmW(:,ind);
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in the SE expressions and
            %update Monte-Carlo estimates 
            signal2_LP_MMSE(:,k,l) = signal2_LP_MMSE(:,k,l) + abs(Hallj_mmW'*w).^2/nbrOfRealizations;            
            signal_LP_MMSE(:,k,l) = signal_LP_MMSE(:,k,l) + Hallj_mmW'*w/nbrOfRealizations;           
            scaling_LP_MMSE(l,k) = scaling_LP_MMSE(l,k) + sum(abs(w).^2,1)/nbrOfRealizations;            
        end
        for ind = 1:length(servedUEs_sub6)
            
            %Extract UE index
            % k = servedUEs(ind);
            k = servedUEs_sub6(ind);
            
            %Normalize LP-MMSE precoding
            w = V_LP_MMSE_sub6(:,ind);
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in the SE expressions and
            %update Monte-Carlo estimates 
            signal2_LP_MMSE(:,k,l) = signal2_LP_MMSE(:,k,l) + abs(Hallj_sub6'*w).^2/nbrOfRealizations;            
            signal_LP_MMSE(:,k,l) = signal_LP_MMSE(:,k,l) + Hallj_sub6'*w/nbrOfRealizations;           
            scaling_LP_MMSE(l,k) = scaling_LP_MMSE(l,k) + sum(abs(w).^2,1)/nbrOfRealizations;            
        end
    end            
end