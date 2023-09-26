function [gainOverNoisedB_col,R_col,pilotIndex_col,D_col,APposition_col,distances_col] = generateSetup_col(L_mmW,L,K_mmW,K,N,coverageRange,coverageRange_sub6,tau_p,nbrOfSetups,seed,APpositions,UEpositions,ASD_varphi,ASD_theta)
%This function generates realizations of the simulation setup described in
%Section 5.3.
%
%INPUT:
%L               = Number of APs per setup
%K               = Number of UEs in the network
%N               = Number of antennas per AP
%tau_p           = Number of orthogonal pilots
%nbrOfSetups     = Number of setups with random UE and AP locations
%seed            = Seed number of pseudorandom number generator
%ASD_varphi      = Angular standard deviation in the local scattering model 
%                  for the azimuth angle (in radians)
%ASD_theta       = Angular standard deviation in the local scattering model
%                  for the elevation angle (in radians)
%
%OUTPUT:
%gainOverNoisedB = Matrix with dimension L x K x nbrOfSetups where
%                  element (l,k,n) is the channel gain (normalized by the
%                  noise variance) between AP l and UE k in setup n
%R               = Matrix with dimension N x N x L x K x nbrOfSetups
%                  where (:,:,l,k,n) is the spatial correlation matrix
%                  between AP l and UE k in setup n, normalized by noise
%pilotIndex      = Matrix with dimension K x nbrOfSetups containing the
%                  pilots assigned to the UEs
%D               = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for cell-free setup
%D_small         = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for small-cell setup
%APpositions     = Vector of length L with the AP locations, where the real
%                  part is the horizontal position and the imaginary part
%                  is the vertical position
%UEpositions     = Vector of length K with UE positions, measured in the
%                  same way as APpositions
%distances       = Matrix with same dimension as gainOverNoisedB containing
%                  the distances in meter between APs and UEs
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


%% Define simulation setup

%Set the seed number if it is specified other than zero
if (nargin>9)&&(seed>0)
    rng(seed)
end

% %Size of the coverage area (as a square with wrap-around)
% squareLength = 1000; %meter

%Communication bandwidth (Hz)
B = 20e6;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters for the model in (5.42)
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading in (5.43)
sigma_sf = 4;

%Decorrelation distance of the shadow fading in (5.43)
decorr = 9;

%Height difference between an AP and a UE (in meters)
distanceVertical = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance


%Prepare to save results
gainOverNoisedB_col = zeros(1,K,nbrOfSetups);
R_col = zeros(N*L,N*L,1,K,nbrOfSetups);
distances_col = zeros(1,K,nbrOfSetups);
pilotIndex_col = zeros(K,nbrOfSetups);
D_col = ones(1,K,nbrOfSetups);
masterAPs = ones(K,1); %the indices of master AP of each UE k 

%% Go through all setups
for n = 1:nbrOfSetups
    APposition_col = mean(APpositions);
    %Compute alternative AP locations by using wrap around
    % wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapHorizontal = repmat([-coverageRange_sub6 0 coverageRange_sub6],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APposition_col,[1 length(wrapLocations)]) + wrapLocations;
    
    %Prepare to store shadowing correlation matrix
    shadowCorrMatrix = sigma_sf^2*ones(K,K);
    shadowAPrealizations = zeros(K,1);    
    
    %Add UEs
    for k = 1:K
        
        %Generate a random UE location in the area
        % UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;
        UEposition = UEpositions(k);  
        %Compute distances assuming that the APs are 10 m above the UEs
        [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEposition,size(APpositionsWrapped))),[],2);
        distances_col(1,k,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
        %If this is not the first UE
        if k-1>0
            
            %Compute distances from the new prospective UE to all other UEs
            shortestDistances = zeros(k-1,1);
            
            for i = 1:k-1
                shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
            end
            
            %Compute conditional mean and standard deviation necessary to
            %obtain the new shadow fading realizations, when the previous
            %UEs' shadow fading realization have already been generated.
            %This computation is based on Theorem 10.2 in "Fundamentals of
            %Statistical Signal Processing: Estimation Theory" by S. Kay
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
            meanvalues = term1*shadowAPrealizations(1:k-1);
            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
            
        else %If this is the first UE
            
            %Add the UE and begin to store shadow fading correlation values
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];
            
        end
        
        %Generate the shadow fading realizations
        shadowing = meanvalues + stdvalue*randn;
        
        %Compute the channel gain divided by noise power
        gainOverNoisedB_col(1,k,n) = constantTerm - alpha*log10(distances_col(1,k,n)) + shadowing' - noiseVariancedBm;
        
        
        
        %Update shadowing correlation matrix and store realizations
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
        shadowAPrealizations(k) = shadowing;
        
        %Store the UE position
%         UEpositions(k) = UEposition;
        
        
        %Determine the master AP for UE k by looking for AP with best
        %channel condition
        masterAPs(k) = 1;
        %Assign orthogonal pilots to the first tau_p UEs according to
        %Algorithm 4.1
        if k <= tau_p
            
            pilotIndex_col(k,n) = k;
            
        else %Assign pilot for remaining UEs
            
            %Compute received power to the master AP from each pilot
            %according to Algorithm 4.1
            pilotinterference = zeros(tau_p,1);
            
            for t = 1:tau_p
                
                pilotinterference(t) = sum(db2pow(gainOverNoisedB_col(pilotIndex_col(1:k-1,n)==t,n)));
                
            end
            
            %Find the pilot with the least receiver power according to
            %Algorithm 4.1
            [~,bestpilot] = min(pilotinterference);
            pilotIndex_col(k,n) = bestpilot;
            
        end

        %Compute nominal angle between UE k and AP 1
        angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(whichpos)); %azimuth angle
        angletoUE_theta = asin(distanceVertical/distances_col(1,k,n));  %elevation angle
        %Generate spatial correlation matrix using the local
        %scattering model in (2.18) and Gaussian angular distribution
        %by scaling the normalized matrices with the channel gain
        if nargin>12
            R_col(:,:,1,k,n) = db2pow(gainOverNoisedB_col(1,k,n))*functionRlocalscattering_mod(N*L,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
        else
            R_col(:,:,1,k,n) = db2pow(gainOverNoisedB_col(1,k,n))*eye(L*N);  %If angular standard deviations are not specified, set i.i.d. fading
        end        
    end        
end
end