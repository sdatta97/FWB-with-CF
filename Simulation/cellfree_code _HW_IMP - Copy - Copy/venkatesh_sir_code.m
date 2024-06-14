function [channelGain_over_noise,R,g_LOS,K_Rician,probLOS] = channel_cellfree_GUE_Krician_non_zero_mean(UE,AP,N,ASD_VALUE,ASD_CORR,rayleigh,BETA,K_factor,cov_area)
%% Noise and channel modelling Constants

% cell-free vs small-cell paper  % eq (52)
D             = cov_area; % 1000/1000; % Side of the square coverage area in Km
h_UE_GUE      = 1.5; % in m   %GUE
h_AP          = 12.5; % in m
f_c           = 1.9*10^3; %Carrier frequency in MHz
sigma_sf_LOS  = 4;
sigma_sf_NLOS = 10; %10;

Band           = 20e6; % Bandwidth
Noise_fig      = 7; % in dB
Noise_temp     = 290; % in kelvin (K)
boltz_const    = 1.381*1e-23; % boltzmann constant
delta          = 0.5; % between 0 and 1. for correlated shadow fading 0 => only at the UE side and 1 => only at the AP side.
antennaSpacing =  1/2; %Antenna spacing
d_DECORR       = 0.1; %KM;   %100; %M  %shadow fading correlation
decorr         = d_DECORR;
%% Noise model
Noise_var     = Band * boltz_const * Noise_temp * Noise_fig;
%Noise_var_dBm = 10*log10(Noise_var) + 30;
Noise_var_dBm = -93.9; %dBm % Noise_fig =7;
%Noise_var_dBm = -94;
%% Simulation area model
AP_locations = (rand(AP,1)*D) + 1i*(rand(AP,1)*D);% randomly generating AP locations within DxD square area.
UE_locations = (rand(UE,1)*D) + 1i*(rand(UE,1)*D);% randomly generating AP locations within DxD square area.  % UAVs, GUEs

wrap_X = repmat([-D 0 D],[3 1]); % wrapping matrix used to wrap AP locations in X direction
wrap_Y = wrap_X';% wrapping matrix used to wrap AP locations in Y direction

wrap_locations       = wrap_X(:)' + 1i*wrap_Y(:)';
AP_locations_wrapped = repmat(AP_locations,[1 length(wrap_locations)]) + repmat(wrap_locations,[AP 1]); % each row of matrix corresponds to AP and its 8 neighbours locations
UE_locations_wrapped = repmat(UE_locations,[1 length(wrap_locations)]) + repmat(wrap_locations,[UE 1]);

%Compute the correlation matrices for the shadow fading   % L=AP, K=UE
shadow_corr_matrix_APs = zeros(AP,AP);
shadow_corr_matrix_UEs = zeros(UE,UE);

for l = 1:AP
    distancetoAP                = min(abs(AP_locations_wrapped - repmat(AP_locations(l),size(AP_locations_wrapped))),[],2);   %2nd dimnsion
    shadow_corr_matrix_APs(:,l) = 2.^(-distancetoAP/decorr);
end
for k = 1:UE
    distancetoUE                = min(abs(UE_locations_wrapped  - repmat(UE_locations(k),size(UE_locations_wrapped ))),[],2);
    shadow_corr_matrix_UEs(:,k) = 2.^(-distancetoUE/decorr);
end
%%
% Prepare to calculate the distance and angles between each UE and AP
UE_AP_dist  = zeros(UE,AP);
UE_AP_angle = zeros(UE,AP);
UE_AP_3D    = zeros(UE,AP);
PL_dash     = zeros(UE,AP);

for ue = 1:UE
    for ap = 1:AP
        [UE_AP_dist(ue,ap), whichAP] = min( abs( UE_locations(ue)*ones(1,length(wrap_locations))-AP_locations_wrapped(ap,:) ) ); %taking the minimum distance out of the 9 virtual AP locations
        %cal 3D dist
        UE_AP_3D(ue,ap)    = sqrt(sum(([real(UE_locations(ue)),imag(UE_locations(ue)),h_UE_GUE/1000]-[real(AP_locations_wrapped(ap,whichAP)), imag(AP_locations_wrapped(ap,whichAP)),h_AP/1000]).^2));
        PL_dash(ue,ap)     = 32.4 + 20*log10(f_c) + 20*log10(UE_AP_3D(ue,ap));  % distance--KM, f_c--MHz
        UE_AP_angle(ue,ap) = angle(UE_locations(ue)-AP_locations_wrapped(ap,whichAP));
    end
end
%% LSF Pathloss model

%Calculate shaddow fading
AP_AP_dist = zeros(AP,AP);
UE_UE_dist = zeros(UE,UE);
for ap1 = 1:AP
    for ap2 = 1:AP
        AP_AP_dist(ap1,ap2) = min(abs(AP_locations(ap1)*ones(1,length(wrap_locations))-AP_locations_wrapped(ap2,:))); %taking the minimum distance out of the 9 virtual AP locations;
    end
end
for ue1 = 1:UE
    for ue2 = 1:UE
        UE_UE_dist(ue1,ue2) = min(abs(UE_locations(ue1)*ones(1,length(wrap_locations))-UE_locations_wrapped(ue2,:))); %taking the minimum distance out of the 9 virtual UE locations
    end
end
Cov_A       = 2.^(-1*AP_AP_dist./d_DECORR);
Cov_B       = 2.^(-1*UE_UE_dist./d_DECORR);
rand_AP     = (randn(AP,1));
rand_UE     = (randn(UE,1));
shadowing_a = sqrtm(Cov_A)*rand_AP;
shadowing_b = sqrtm(Cov_B)*rand_UE;

%% ------------------------------------------
R_norm        = zeros(N,N,UE,AP,length(ASD_VALUE));
maxdistLOS    = 300;
probLOS       = zeros(UE,AP);
shaddow_fad   = zeros(UE,AP);
path_loss_dBm = zeros(UE,AP);
%---------------------------------------
% K_Rician = zeros(UE,AP);
%---------------------------------------
K_Rician=10.^((13-0.03*UE_AP_dist*1000)./10);
%K_Rician=20*ones(UE,AP);
%---------------------------------------
%Calculate Pathloss component with shaddow fading, channel mean and covariance matrices
for ue = 1:UE
    for ap = 1:AP
        %-------------------------------------------------------
        if rayleigh==1  %rayleigh
            probLOS(ue,ap) = 0; K_Rician(ue,ap) = 0;
            shaddow_fad(ue,ap)  = sigma_sf_NLOS*(sqrt(delta)*shadowing_a(ap) + sqrt(1-delta)*shadowing_b(ue));% correlated shadow fading model from paper
            path_loss_dBm(ue,ap) = -34.53-38*log10(UE_AP_dist(ue,ap)*1000) + shaddow_fad(ue,ap);  %Pathloss component for Rayleigh fading
        else
            %     Generate shadow fading realizations
            %-------------------------------------------------------
            %             if (UE_AP_dist(ue,ap)*1000) <= 300
            %                 probLOS(ue,ap)      = ((maxdistLOS-(UE_AP_dist(ue,ap)*1000))./maxdistLOS);
            %                 K_Rician(ue,ap)     = 10.^((13-0.03*UE_AP_dist(ue,ap)*1000)./10);
            %                 shaddow_fad(ue,ap)  = sigma_sf_LOS*(sqrt(delta)*shadowing_a(ap)+ sqrt(1-delta)*shadowing_b(ue));% correlated shadow fading model from paper
            %                 path_loss_dBm(ue,ap) = -30.18-26*log10(UE_AP_dist(ue,ap)*1000) +shaddow_fad(ue,ap);  %Pathloss component for Rician Fading
            %             else
            %                 probLOS(ue,ap) = 0; K_Rician(ue,ap) = 0;
            %                 shaddow_fad(ue,ap)   = sigma_sf_NLOS*(sqrt(delta)*shadowing_a(ap)+ sqrt(1-delta)*shadowing_b(ue));% correlated shadow fading model from paper
            %                 path_loss_dBm(ue,ap) = -34.53-38*log10(UE_AP_dist(ue,ap)*1000) +shaddow_fad(ue,ap);  %Pathloss component for Rayleigh fading
            %             end
            %-------------------------------------------------------
            if K_Rician(ue,ap)~=0
                probLOS(ue,ap)=((maxdistLOS-(UE_AP_dist(ue,ap)*1000))./maxdistLOS);
                shaddow_fad(ue,ap)    = sigma_sf_LOS*(sqrt(delta)*shadowing_a(ap)+ sqrt(1-delta)*shadowing_b(ue));% correlated shadow fading model from paper
                path_loss_dBm(ue,ap) = -30.18-26*log10(UE_AP_dist(ue,ap)*1000)+shaddow_fad(ue,ap);  %Pathloss component for Rician Fading
            elseif K_Rician(ue,ap)==0
                probLOS(ue,ap)= 0;
                shaddow_fad(ue,ap)    = sigma_sf_NLOS*(sqrt(delta)*shadowing_a(ap)+ sqrt(1-delta)*shadowing_b(ue));% correlated shadow fading model from paper
                path_loss_dBm(ue,ap) = -34.53-38*log10(UE_AP_dist(ue,ap)*1000)+shaddow_fad(ue,ap);  %Pathloss component for Rayleigh fading
            end
            %-------------------------------------------------------
        end
        %---------------CORR
        if ASD_CORR ==1
            for iASD = 1:length(ASD_VALUE)
                if ASD_VALUE(iASD) ~=0
                    R_norm(:,:,ue,ap,iASD) = functionRlocalscattering(N,UE_AP_angle(ue,ap),ASD_VALUE(iASD),0.5,'Gaussian'); %Normalized covariance matrix (Correlated)
                elseif ASD_VALUE(iASD) == 0
                    R_norm(:,:,ue,ap,iASD) = eye(N);                                                                 %Normalized covariance matrix (Uncorrelated)
                end
            end
        end
    end
end

if ASD_CORR ==0   % EXP-CORR
    R_EXP = zeros(N,N,length(ASD_VALUE));
    for iASD1 = 1:length(ASD_VALUE)
        for u1 = 1:N
            for u2 = 1:N
                R_EXP(u1,u2,iASD1) = ASD_VALUE(iASD1)^(abs(u1-u2));
            end
        end
    end
end
%% Channel gain over noise
if BETA == 0
    channelGain_over_noise = 10.^(0.1*(path_loss_dBm-Noise_var_dBm)); %10*rand(UE,AP); % ones(UE,AP);    (Linear scale)

elseif BETA==1
    channelGain_over_noise = 10*ones(UE,AP); % 10*rand(UE,AP); %ones(UE,AP);
    probLOS                = ones(UE,AP); %randi(2,UE,AP)-1;
    K_Rician               = db2pow(K_factor)*ones(UE,AP).*probLOS;
end

%% Calculate --- Covariance matrices
R = zeros(N,N,AP,UE,length(ASD_VALUE));
for ap=1:AP
    for ue=1:UE
        %-----------CORR
        for iASD = 1:length(ASD_VALUE)
            if K_Rician(ue,ap) == 0
                if ASD_CORR == 1
                    R(:,:,ap,ue,iASD) = channelGain_over_noise(ue,ap)*R_norm(:,:,ue,ap,iASD);
                elseif ASD_CORR == 0
                    R(:,:,ap,ue,iASD) = channelGain_over_noise(ue,ap)*R_EXP(:,:,iASD);
                end
            else
                if ASD_CORR == 1
                    R(:,:,ap,ue,iASD) = (1/(K_Rician(ue,ap)+1))*channelGain_over_noise(ue,ap)*R_norm(:,:,ue,ap,iASD);
                elseif ASD_CORR == 0
                    R(:,:,ap,ue,iASD) = (1/(K_Rician(ue,ap)+1))*channelGain_over_noise(ue,ap)*R_EXP(:,:,iASD);
                end
            end
        end
    end
end
%% LOS component
g_LOS = zeros(N,AP,UE);
for l = 1:AP
    for k = 1:UE
        if K_Rician(k,l)== 0
            g_LOS(:,l,k)  = zeros(N,1);
        else
            g_LOS(:,l,k) = sqrt(channelGain_over_noise(k,l)/(K_Rician(k,l)+1) )* sqrt(K_Rician(k,l))*(exp(1i*2*pi.*(0:(N-1))*sin(UE_AP_angle(k,l))*antennaSpacing));      %Normalized Mean vector
        end
    end
end
%%
% scatter(real(AP_locations), imag(AP_locations), 'o');
% hold on;
% scatter(real(UE_locations), imag(UE_locations), '+');
% hold off;
% legend('AP','UE')
% title('Simulation area layout');
% xlabel('D [Km]');
% ylabel('D [Km]');
end

