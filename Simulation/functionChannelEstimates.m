function [Hhat_mmW,Hhat_sub6,H_mmW,H_sub6,B,C] = functionChannelEstimates(R,nbrOfRealizations,L_mmW,L,K_mmW,K,N,N_UE_mmW,tau_p,pilotIndex,p)
%Generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are assumed to be correlated
%Rayleigh fading and the MMSE estimator is used.
%
%INPUT:
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k,
%                    normalized by noise variance
%nbrOfRealizations = Number of channel realizations
%L                 = Number of APs
%K                 = Number of UEs in the network
%N                 = Number of antennas per AP
%tau_p             = Number of orthogonal pilots
%pilotIndex        = Vector containing the pilot assigned to each UE
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%Hhat         = Matrix with dimension L*N x nbrOfRealizations x K where
%               (:,n,k) is the estimated collective channel to UE k in
%               channel realization n.
%H            = Matrix with dimension L*N x nbrOfRealizations x K with the
%               true channel realizations. The matrix is organized in the
%               same way as Hhat.
%B            = Matrix with dimension N x N x L x K where (:,:,l,k) is the
%               spatial correlation matrix of the channel estimate between
%               AP l and UE k, normalized by noise variance
%C            = Matrix with dimension N x N x L x K where (:,:,l,k) is the
%               spatial correlation matrix of the channel estimation error
%               between AP l and UE k, normalized by noise variance
%
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




%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
% H = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K));
H_mmW = (randn(L*N,nbrOfRealizations,N_UE_mmW,K_mmW)+1i*randn(L*N,nbrOfRealizations,N_UE_mmW,K_mmW));
H_sub6 = (randn(L*N,nbrOfRealizations,K-K_mmW)+1i*randn(L*N,nbrOfRealizations,K-K_mmW));


%Go through all channels and apply the spatial correlation matrices
for n = 1:nbrOfRealizations
    for l = 1:L
        for k = 1:K_mmW            
            %Apply correlation to the uncorrelated channel realizations
            Rsqrt = sqrtm(R(:,:,l,k));
            % H_mmW((l-1)*N+1:l*N,n,:,k) = sqrt(0.5)*Rsqrt*H_mmW((l-1)*N+1:l*N,n,:,k);
            H_mmW((l-1)*N+1:l*N,n,:,k) = sqrt(0.5)*Rsqrt*reshape(H_mmW((l-1)*N+1:l*N,n,:,k),[N,N_UE_mmW]);           
        end
        for k = 1:(K-K_mmW)        
            %Apply correlation to the uncorrelated channel realizations
            Rsqrt = sqrtm(R(:,:,l,k+K_mmW));
            % H_sub6((l-1)*N+1:l*N,n,:,k-K_mmW) = sqrt(0.5)*Rsqrt*H_sub6((l-1)*N+1:l*N,n,:,k-K_mmW);
            H_sub6((l-1)*N+1:l*N,n,k) = sqrt(0.5)*Rsqrt*H_sub6((l-1)*N+1:l*N,n,k);            
        end
    end
end


%% Perform channel estimation

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));

%Prepare to store results
% Hhat = zeros(L*N,nbrOfRealizations,K);
Hhat_mmW = zeros(L*N,nbrOfRealizations,N_UE_mmW,K_mmW);
Hhat_sub6 = zeros(L*N,nbrOfRealizations,K-K_mmW);

if nargout>2
    B = zeros(size(R));
end

if nargout>3
    C = zeros(size(R));
end


%Go through all APs
for l = 1:L
    
    %Go through all pilots
    for t = 1:tau_p
        
        %Compute processed pilot signal for all UEs that use pilot t
        %according to (4.4) with an additional scaling factor \sqrt{tau_p}
        % yp = sqrt(p)*tau_p*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + sqrt(tau_p)*Np(:,:,l,t);
        bool_pilots = (t==pilotIndex);
        yp_mmW = sqrt(p)*tau_p*sum(H_mmW((l-1)*N+1:l*N,:,:,bool_pilots(1:K_mmW)),4) + sqrt(tau_p)*Np(:,:,l,t);
        yp_sub6 = sqrt(p)*tau_p*sum(H_sub6((l-1)*N+1:l*N,:,bool_pilots(1+K_mmW:K)),3) + sqrt(tau_p)*Np(:,:,l,t);
        
        %Compute the matrix in (4.6) that is inverted in the MMSE estimator
        %in (4.5)
        PsiInv = (p*tau_p*sum(R(:,:,l,t==pilotIndex),4) + eyeN);
        %Go through all UEs that use pilot t
        for k = find(t==pilotIndex)'           
            %Compute the MMSE estimate
            RPsi = R(:,:,l,k) / PsiInv;
            % Hhat((l-1)*N+1:l*N,:,k) = sqrt(p)*RPsi*yp;
            for n=1:nbrOfRealizations
                if (k<=K_mmW)
                    Hhat_mmW((l-1)*N+1:l*N,n,:,k) = sqrt(p)*RPsi*reshape(yp_mmW(:,n,:),[N,N_UE_mmW]);
                else
                    Hhat_sub6((l-1)*N+1:l*N,n,k-K_mmW) = sqrt(p)*RPsi*yp_sub6(:,n);
                end
            end
            if nargout>2
                %Compute the spatial correlation matrix of the estimate
                %according to (4.7)
                B(:,:,l,k) = p*tau_p*RPsi*R(:,:,l,k);
            end
            
            if nargout>3
                %Compute the spatial correlation matrix of the estimation
                %error according to (4.9)
                C(:,:,l,k) = R(:,:,l,k) - B(:,:,l,k);
            end
            
        end
    end    
end
