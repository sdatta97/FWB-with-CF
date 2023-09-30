function [SE_eq, SE] = functionDownlinkSE_sumSE_distv3_eq(beta,preLogFactor,L,K,K_mmW,N_AP,N_UE,D,rhomax,tau_p)
%Compute downlink SE according to Corollary 6.3 with sum SE maximizing power allocation
%in Algorithm 7.6
%
%INPUT:
%bk             = Matrix with dimension L x K where (1:La,k) is the non-zero
%                 portion of \tilde{b}_k in (7.25) (La is the number of
%                 APs serving UE k)
%Ck             = Matrix with dimension L x L x K x K where (1:La,1:La,k,i) is the 
%                 non-zero portion of \tilde{C}_{ki} in (7.26) (La is the
%                 number of APs serving UE i)
%preLogFactor   = Pre-log factor in the SE expression of Corollary 6.3
%L              = Number of APs
%K              = Number of UEs in the network
%rhomax         = Maximum allowable AP transmit power
%tau_p          = Length of pilot sequences
%
%OUTPUT:
%SE             = SEs achieved with sum SE maximizing power allocation algorithm
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

%Initialize the current objective function
objec_new = 0; % inf;

%Prepare array to store the number of APs serving a specficic UE
La = zeros(K,1);
%Prepare cell to store the AP indices serving a specficic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ = cell(K,1);
p_fac= 1; %ratio of mmW to sub-6 powers
%Construc the above array and cells
for k = 1:K
    servingAPs = find(D(:,k)==1);
    NoservingAPs = find(D(:,k)==0);
    
    Serv{k} = servingAPs;
    NoServ{k} = NoservingAPs;
    
    La(k) = length(servingAPs);
    beta(:,k) = beta(:,k).*D(:,k);
end
beta_opt = zeros(sum(La),K);
for k=1:K
    if (k==1)
        beta_opt(1:La(k),k) = beta(Serv{k},k);
    else
        beta_opt(1+sum(La(1:k-1)):sum(La(1:k)),k) = beta(Serv{k},k);
    end
end
%Intialize the difference between current and previous objective values 
diff = 100;
%Initialize iterates
eta_eq = zeros(L,K);
if (K_mmW == 0)
    for l = 1:L
        for k = 1:K
            if ismember(l,Serv{k})
                eta_eq(l,k) = 1./(N_AP*N_UE*sum(beta(l,:)));
            end
        end
    end
else
    for l = 1:L
        for k = 1:K
            if ismember(l,Serv{k})
                if (k<=K_mmW)
                    eta_eq(l,k) = p_fac./(N_AP*N_UE*(p_fac*beta(l,1)+sum(beta(l,2:K))));
                else
                    eta_eq(l,k) = 1./(N_AP*N_UE*(p_fac*beta(l,1)+sum(beta(l,2:K))));
                end
            end
        end
    end
end
lambda_eq = zeros(K,1); %sum((sqrt(eta_eq)*D).*beta,1)';
zeta_eq = zeros(K,1);
for k = 1:K
    lambda_eq(k) = (sqrt(eta_eq(:,k)))'*beta(:,k);
    zeta_eq(k) = (lambda_eq(k)^2)/(1/(rhomax*N_AP*N_AP) + (N_UE/N_AP)*beta(:,k)'*sum(beta.*eta_eq,2));
end
lambda_old = lambda_eq;
zeta_old = zeta_eq;
% eta = eta_eq;
SE_eq = preLogFactor*log(1+zeta_eq)/log(2);
SE=SE_eq;
end
