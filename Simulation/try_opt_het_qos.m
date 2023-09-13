%This Matlab script can be used to reproduce Figure 7.3 in the monograph:
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

%Empty workspace and close figures
close all;
clear;


%% Define simulation setup
lambda_BS = 10; %25;
lambda_UE = 0; %10;
lambda_UE_sub6 = 10; %50; %:10:50;
coverageRange = 100;
length_area = 2*coverageRange;   
width_area = 2*coverageRange;
coverageRange_sub6 = 1000;
length_area_sub6 = 2*coverageRange_sub6;   
width_area_sub6 = 2*coverageRange_sub6;
%height receiver (UE), approximately the height a human holds the phone
hr = 1.4;
ht = 5;
ht_sub6 = 4;
areaDimensions = [width_area, length_area, ht];
areaDimensions_sub6 = [width_area_sub6, length_area_sub6, ht_sub6];


%Number of setups with random UE locations
nbrOfSetups = 10;


        
      
%Number of channel realizations per setup
nbrOfRealizations = 10;

% %Number of APs in the cell-free network
% L = 10;

%Number of antennas per AP
N = 10;

%Number of antennas per UE
N_UE_mmW = 8;
N_UE_sub6 = 2;

%Number of UEs in the network
% K = 40;

%Length of the coherence block
tau_c = 200;

%Compute number of pilots per coherence block
tau_p = 10;

%Compute the prelog factor assuming only downlink data transmission
preLogFactor = (tau_c-tau_p)/tau_c;

%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15); %azimuth angle
ASD_theta = deg2rad(15);  %elevation angle

%Total uplink transmit power per UE (mW)
p = 100;

%Total downlink transmit power per AP (mW)
rho_tot = 200;

%min-QoS reqs
rmin = 1e9;
rmin_sub6 = 1e7;
%Prepare to save simulation results

% SE_DL_LPMMSE_equal = zeros(K,nbrOfSetups); %Equal
% SE_DL_LPMMSE_fractional = zeros(K,nbrOfSetups); %FPA, \upsilon = 0.5
% SE_DL_LPMMSE_fractional2 = zeros(K,nbrOfSetups); %FPA, \upsion = -0.5
% SE_DL_LPMMSE_maxmin = zeros(K,nbrOfSetups); %MMF
% SE_DL_LPMMSE_sumSE = zeros(K,nbrOfSetups); %SumSE
SE_DL_LPMMSE_equal_mean = zeros(nbrOfSetups); %Equal
SE_DL_LPMMSE_fractional_mean = zeros(nbrOfSetups); %FPA, \upsilon = 0.5
SE_DL_LPMMSE_maxmin_mean = zeros(nbrOfSetups); %MMF
SE_DL_LPMMSE_sumSE_mean = zeros(nbrOfSetups); %SumSE

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    nn = poissrnd(lambda_BS*pi*(coverageRange/1000)^2);
    while (nn==0)
        nn = poissrnd(lambda_BS*pi*(coverageRange/1000)^2);       
    end
    L_mmW = nn;
    L = poissrnd(lambda_BS*pi*(coverageRange_sub6/1000)^2);
    %%UE locations
    %%UE location
    % n = poissrnd(lambda_UE*pi*(params.coverageRange_sub6/1000)^2);
    % while (n==0)
    %     n = poissrnd(lambda_UE*pi*(params.coverageRange_sub6/1000)^2);
    % end
    % params.numUE = n;
    K_mmW = 1;
    K_sub6 = poissrnd(lambda_UE_sub6*pi*(coverageRange_sub6/1000)^2);
    K = K_mmW+K_sub6;
    SE_DL_LPMMSE_equal = zeros(K,1); %Equal
    SE_DL_LPMMSE_fractional = zeros(K,1); %FPA, \upsilon = 0.5
    SE_DL_LPMMSE_maxmin = zeros(K,1); %MMF
    SE_DL_LPMMSE_sumSE = zeros(K,1); %SumSE
    %Generate one setup with UEs at random locations
    [gainOverNoisedB,R,pilotIndex,D,D_small] = generateSetup(L_mmW,L,K_mmW,K,N,coverageRange,coverageRange_sub6,tau_p,1,0,ASD_varphi,ASD_theta);
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cell-free APs
    [Hhat_mmW,Hhat_sub6,H_mmW,H_sub6,B,C] = functionChannelEstimates(R,nbrOfRealizations,L_mmW,L,K_mmW,K,N,N_UE_mmW,N_UE_sub6,tau_p,pilotIndex,p);
    % [~,l_idx] = max(mean(chgains,1));
    % Full uplink power for the computation of precoding vectors using
    % virtual uplink-downlink duality
    p_full = p*ones(K,1);
   
    gainOverNoise = db2pow(gainOverNoisedB);

    %Equal power allocation
    rho_dist_equal = (rho_tot/tau_p)*ones(L,K);

    %Compute the power allocation in (7.47) for distributed precoding
    rho_dist = zeros(L,K); % with exponent 0.5
    
    for l = 1:L
        
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute denominator in (7.47)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));

        for ind = 1:length(servedUEs)
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
        end
        
    end
    
    %Obtain the expectations for the computation of the terms in
    %(7.25)-(7.26)
    %before offload sub-6
    % [signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
    %  functionComputeExpectations(Hhat_sub6,H_sub6,D(:,(1+K_mmW):end),C(:,:,:,(1+K_mmW):end),nbrOfRealizations,N,K-K_mmW,L,p_full((1+K_mmW):end));
     
    %Prepare arrays to store the vectors \tilde{b}_k in (7.25) and matrices
    %\tilde{C}_{ki} in (7.26)
 %    [bk, Ck] = ...
 % functionComputeExpectationsv2(Hhat_mmW, H_mmW, Hhat_sub6, H_sub6,D,C,nbrOfRealizations,N,N_UE_mmW, N_UE_sub6,K,K_mmW,L,L_mmW,p_full);
% [bk_mmW, Ck_mmW, bk_sub6, Ck_sub6] = functionComputeExpectationsv3(Hhat_mmW, H_mmW, Hhat_sub6, H_sub6,D,C,nbrOfRealizations,N,N_UE_mmW, N_UE_sub6,K,K_mmW,L,L_mmW,p_full, gainOverNoise);
[bk, Ck] = functionComputeExpectationsv4(Hhat_mmW, H_mmW, Hhat_sub6, H_sub6,D,C,nbrOfRealizations,N,N_UE_mmW, N_UE_sub6,K,K_mmW,L,L_mmW,p_full, gainOverNoise);

    % %Go through all UEs
    % for k = 1:K
    %     %Find the APs that serve UE k
    %     servingAPs = find(D(:,k)==1);
    %     %The number of APs that serve UE k
    %     La = length(servingAPs);
    %     %Compute the vector in (7.25) for UE k (only the non-zero indices correspondig to 
    %     %serving APs are considered)
    %     bk(1:La,k) = real(vec(signal_LP_MMSE(k,k,servingAPs)))./sqrt(scaling_LP_MMSE(servingAPs,k));
    % 
    %     %Go through all UEs
    %     for i = 1:K
    %         %Find the APs that serve UE i
    %         servingAPs = find(D(:,i)==1);
    %         %The number of APs that serve UE i
    %         La = length(servingAPs);
    %         %Compute the matrices in (7.26) (only the non-zero indices are
    %         %considered)
    %         if i==k
    %            Ck(1:La,1:La,k,k) = bk(1:La,k)*bk(1:La,k)';
    %         else
    %            Ck(1:La,1:La,k,i) = diag(1./sqrt(scaling_LP_MMSE(servingAPs,i)))...
    %                *(vec(signal_LP_MMSE(k,i,servingAPs))...
    %                *vec(signal_LP_MMSE(k,i,servingAPs))')...
    %                *diag(1./sqrt(scaling_LP_MMSE(servingAPs,i)));
    %         end            
    %         for j = 1:La
    %             Ck(j,j,k,i) = signal2_LP_MMSE(k,i,servingAPs(j))/scaling_LP_MMSE(servingAPs(j),i);
    %         end
    %     end
    % end
    
    %Take the real part (in the SINR expression,the imaginary terms cancel
    %each other)
    bk = real(bk);
    Ck = real(Ck);
%     bk_mmW = real(bk_mmW);
%     Ck_mmW = real(Ck_mmW);
%     bk_sub6 = real(bk_sub6);
%     Ck_sub6 = real(Ck_sub6);
    %Compute hte square roots of the power allocation coefficients
    %corresponding to (7.24)
    tilrho = sqrt(rho_dist_equal);
    tilrho1 = sqrt(rho_dist);

    %Go through all UEs
    for k = 1:K-K_mmW
        %Find APs that serve UE k
        servingAPs = find(D(:,k+K_mmW)==1);
        %The number of APs that serve UE k
        La = length(servingAPs);
        for nn = 1:N_UE_sub6
            %Compute the numerator and denominator of (7.23) for equal and FPA
            %schemes with two different exponents
            numm = abs(bk(1:La,k+K_mmW)'*tilrho(servingAPs,k+K_mmW))^2;
%             numm = abs(reshape(bk_sub6(servingAPs,k,nn),[La,1])'*tilrho(servingAPs,k+K_mmW))^2;
            denomm = 1-numm;
            
            numm1 = abs(bk(1:La,k+K_mmW)'*tilrho1(servingAPs,k+K_mmW))^2;
%             numm1 = abs(reshape(bk_sub6(servingAPs,k,nn),[La,1])'*tilrho1(servingAPs,k+K_mmW))^2;
            denomm1 = 1-numm1;

            for i = 1:K-K_mmW
                servingAPs = find(D(:,i+K_mmW)==1);
                La = length(servingAPs);
                denomm = denomm+tilrho(servingAPs,i+K_mmW)'*Ck(1:La,1:La,k+K_mmW,i+K_mmW)*tilrho(servingAPs,i+K_mmW);
                denomm1 = denomm1+tilrho1(servingAPs,i+K_mmW)'*Ck(1:La,1:La,k+K_mmW,i+K_mmW)*tilrho1(servingAPs,i+K_mmW);        
%                 denomm = denomm+tilrho(servingAPs,i+K_mmW)'*reshape(Ck_sub6(servingAPs,servingAPs,k,i,nn),[La,La])*tilrho(servingAPs,i+K_mmW);
%                 denomm1 = denomm1+tilrho1(servingAPs,i+K_mmW)'*reshape(Ck_sub6(servingAPs,servingAPs,k,i,nn),[La,La])*tilrho1(servingAPs,i+K_mmW);        
            end
            %Compute SEs using SINRs in (7.23) and Corollary 6.3 for equal and
            %FPA schemes with two different exponents
            SE_DL_LPMMSE_equal(k+K_mmW,n) = SE_DL_LPMMSE_equal(k+K_mmW,n) + preLogFactor*log2(1+numm/denomm);
            SE_DL_LPMMSE_fractional(k+K_mmW,n) = SE_DL_LPMMSE_fractional(k+K_mmW,n) + preLogFactor*log2(1+numm1/denomm1);
        end
    end
    
    %Compute SE according to Corollary 6.3 with max-min fair power
    %allocation in Algorithm 7.5
%     SE_DL_LPMMSE_maxmin((1+K_mmW):end,n) = functionDownlinkSE_maxmin_dist(bk_sub6,Ck_sub6,preLogFactor,L,K-K_mmW,D(:,(1+K_mmW):end),rho_tot);  
    %Compute SE according to Corollary 6.3 with sum SE maximizing power
    %allocation in Algorithm 7.6
    SE_DL_LPMMSE_sumSE((1+K_mmW):end,n) =  functionDownlinkSE_sumSE_dist(bk(:,(1+K_mmW):end),Ck(:,:,(1+K_mmW):end,(1+K_mmW):end),preLogFactor,L,K-K_mmW,D(:,(1+K_mmW):end),rho_tot,tau_p);   
%     SE_DL_LPMMSE_sumSE((1+K_mmW):end,n) =  functionDownlinkSE_sumSE_dist(bk_sub6,Ck_sub6,preLogFactor,L,K-K_mmW,D(:,(1+K_mmW):end),rho_tot,tau_p);   
% SE_DL_LPMMSE_sumSE((1+K_mmW):end,n) =  functionDownlinkSE_sumSE_distv2(bk_sub6,Ck_sub6,preLogFactor,L,K-K_mmW,N_UE_sub6,D(:,(1+K_mmW):end),rho_tot,tau_p);   

    %% 
    %excluding mmW serving gNB
    % Full uplink power for the computation of precoding vectors using
    % virtual uplink-downlink duality
    D(l_idx,(1+K_mmW):end) = 0;
    [Hhat_mmW,Hhat_sub6,H_mmW,H_sub6,B,C] = functionChannelEstimates(R,nbrOfRealizations,L_mmW,L,K_mmW,K,N,N_UE_mmW,tau_p,pilotIndex,p);
    p_full = p*ones(K,1);
   
    gainOverNoise = db2pow(gainOverNoisedB);

    %Equal power allocation
    rho_dist_equal = (rho_tot/tau_p)*ones(L,K);

    %Compute the power allocation in (7.47) for distributed precoding
    rho_dist = zeros(L,K); % with exponent 0.5
    for l = 1:L
        %Extract which UEs are served by AP l
        servedUEs = find(D(l,:)==1);
        
        %Compute denominator in (7.47)
        normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));

        for ind = 1:length(servedUEs)
            rho_dist(l,servedUEs(ind)) = rho_tot*sqrt(gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
        end        
    end
    
    %Obtain the expectations for the computation of the terms in
    %(7.25)-(7.26)   
    %after offload sub-6
     [signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
     functionComputeExpectations(Hhat_sub6,H_sub6,D(:,(1+K_mmW):end),C(:,:,:,(1+K_mmW):end),nbrOfRealizations,N,K-K_mmW,L,p_full((1+K_mmW):end));

    %Prepare arrays to store the vectors \tilde{b}_k in (7.25) and matrices
    %\tilde{C}_{ki} in (7.26)
    bk = zeros(L,K-K_mmW);
    Ck = zeros(L,L,K-K_mmW,K-K_mmW);
    %Go through all UEs
    for k = 1:K-K_mmW
        %Find the APs that serve UE k
        servingAPs = find(D(:,k+K_mmW)==1);
        %The number of APs that serve UE k
        La = length(servingAPs);
        %Compute the vector in (7.25) for UE k (only the non-zero indices correspondig to 
        %serving APs are considered)
        bk(1:La,k) = real(vec(signal_LP_MMSE(k,k,servingAPs)))./sqrt(scaling_LP_MMSE(servingAPs,k));
        
        %Go through all UEs
        for i = 1:K-K_mmW
            %Find the APs that serve UE i
            servingAPs = find(D(:,i+K_mmW)==1);
            %The number of APs that serve UE i
            La = length(servingAPs);          
           %Compute the matrices in (7.26) (only the non-zero indices are
            %considered)
            if i==k
               Ck(1:La,1:La,k,k) = bk(1:La,k)*bk(1:La,k)';
            else
                   Ck(1:La,1:La,k,i) = diag(1./sqrt(scaling_LP_MMSE(servingAPs,i)))...
                   *(vec(signal_LP_MMSE(k,i,servingAPs))...
                   *vec(signal_LP_MMSE(k,i,servingAPs))')...
                   *diag(1./sqrt(scaling_LP_MMSE(servingAPs,i)));
            end            
            for j = 1:La
                Ck(j,j,k,i) = signal2_LP_MMSE(k,i,servingAPs(j))/scaling_LP_MMSE(servingAPs(j),i);
            end
        end
    end
    
    %Take the real part (in the SINR expression,the imaginary terms cancel
    %each other)
    Ck = real(Ck);
    
    %Compute hte square roots of the power allocation coefficients
    %corresponding to (7.24)
    tilrho = sqrt(rho_dist_equal);
    tilrho1 = sqrt(rho_dist);

    %Go through all UEs
    % for k = 1:K
    for k = 1:K-K_mmW
        %Find APs that serve UE k
        servingAPs = find(D(:,k+K_mmW)==1);
        %The number of APs that serve UE k
        La = length(servingAPs);
        %Compute the numerator and denominator of (7.23) for equal and FPA
        %schemes with two different exponents
        numm = abs(bk(1:La,k)'*tilrho(servingAPs,k+K_mmW))^2;
        denomm = 1-numm;
        
        numm1 = abs(bk(1:La,k)'*tilrho1(servingAPs,k+K_mmW))^2;
        denomm1 = 1-numm1;
        
        for i = 1:K-K_mmW
            servingAPs = find(D(:,i+K_mmW)==1);
            La = length(servingAPs);
            denomm = denomm+tilrho(servingAPs,i+K_mmW)'*Ck(1:La,1:La,k,i)*tilrho(servingAPs,i+K_mmW);
            denomm1 = denomm1+tilrho1(servingAPs,i+K_mmW)'*Ck(1:La,1:La,k,i)*tilrho1(servingAPs,i+K_mmW);        
        end
        %Compute SEs using SINRs in (7.23) and Corollary 6.3 for equal and
        %FPA schemes with two different exponents
        SE_DL_LPMMSE_equal(k+K_mmW,n) = preLogFactor*log2(1+numm/denomm);
        SE_DL_LPMMSE_fractional(k+K_mmW,n) = preLogFactor*log2(1+numm1/denomm1);
    end
    
    %Compute SE according to Corollary 6.3 with max-min fair power
    %allocation in Algorithm 7.5
    SE_DL_LPMMSE_maxmin((1+K_mmW):end,n) = functionDownlinkSE_maxmin_dist(bk,Ck,preLogFactor,L,K-K_mmW,D(:,(1+K_mmW):end),rho_tot);  
    %Compute SE according to Corollary 6.3 with sum SE maximizing power
    %allocation in Algorithm 7.6
    SE_DL_LPMMSE_sumSE((1+K_mmW):end,n) =  functionDownlinkSE_sumSE_dist(bk,Ck,preLogFactor,L,K-K_mmW,D(:,(1+K_mmW):end),rho_tot,tau_p);   
end

% Plot Figure 7.3
figure;
hold on; box on;
set(gca,'fontsize',16);

plot(sort(SE_DL_LPMMSE_equal(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(SE_DL_LPMMSE_fractional(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',2);
plot(sort(SE_DL_LPMMSE_maxmin(:)),linspace(0,1,K*nbrOfSetups),'b-.','LineWidth',2);
plot(sort(SE_DL_LPMMSE_sumSE(:)),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'Equal', 'FPA, $\upsilon=0.5$', 'MMF','SumSE' },'Interpreter','Latex','Location','SouthEast');
xlim([0 12]);