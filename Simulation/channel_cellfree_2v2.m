
% function [channelGain_over_noise, interAPGain_over_noise, interUEGain_over_noise, sigma_sf, delta] = channel_cellfree_2v2(Kd,Ku,M,Band)
% function [channelGain_over_noise] = channel_cellfree_2v2(Kd,Ku,M,Band, AP_locations, UE_locations)
function [channelGain_over_noise, ricianFactor] = channel_cellfree_2v2(Kd,M,Band, AP_locations, UE_locations)
%% Noise and channel modelling Constants
% K = Ku + Kd;
K = Kd;
D     = 1000; % Side of the square coverage area in m
% D0    = 10; % in m
% D1    = 50; % in m
h_UE  = 1.5; % in m
h_AP  = 15; % in m
f_c   = 1.9*10^3; %Carrier frequency in MHz
%Band  = 20e6; % Bandwidth
Noise_fig_dB   = 9; % in dB
Noise_fig = 10^(0.1*Noise_fig_dB);
Noise_temp  = 290; % in kelvin (K)
boltz_const = 1.381*1e-23; % boltzmann constant
sigma_sf = 8;
delta       = 0.5; % between 0 and 1. for correlated shadow fading 0 => only at the UE side and 1 => only at the AP side.
d_DECORR    = 100; % Decorrelation distance for the shadow fading (in m)
%Correlation model parameters
% CORR        = 0; % switch it to 1 for spatially correlated system model.
% N_AP        = 1; %number of antennas at each AP used for spatial correlation model.
% ASD_deg     = 10;
% antenna_spacing = 0.5;
%% Noise model
Noise_var = Band * boltz_const * Noise_temp * Noise_fig; 
Noise_var_dB = 10*log10(Noise_var);
%% Simulation area model
% AP_locations = (rand(M,1)*D) + 1i*(rand(M,1)*D);% randomly generating AP locations within DxD square area. 
% UE_locations = (rand(K,1)*D) + 1i*(rand(K,1)*D);% randomly generating UE locations within DxD square area. 

% wrap_X = repmat([-D 0 D],[3 1]); % wrapping matrix used to wrap AP locations in X direction
% wrap_Y = wrap_X';% wrapping matrix used to wrap AP locations in Y direction
% 
% wrap_locations       = wrap_X(:)' + 1i*wrap_Y(:)';
% AP_locations_wrapped = repmat(AP_locations,[1 length(wrap_locations)]) + repmat(wrap_locations,[M 1]); % each row of matrix corresponds to AP and its 8 neighbours locations
% UE_locations_wrapped = repmat(UE_locations,[1 length(wrap_locations)]) + repmat(wrap_locations,[K 1]); % each row of matrix corresponds to AP and its 8 neighbours locations

% Prepare to calculate the distance and angles between each UE and AP
UE_AP_dist  = zeros(M,K);
AP_AP_dist  = zeros(M,M);
UE_UE_dist  = zeros(K,K);
% constantTerm  = 46.3+33.9*log10(f_c)-13.82*log10(h_AP)-(1.1*log10(f_c)-0.7)*h_UE+(1.56*log10(f_c)-0.8);
channelGaindB = zeros(M,K);
ricianFactor = zeros(M,K);
for ue = 1:K
    for ap = 1:M
        UE_AP_dist(ap,ue) = sqrt(sum((UE_locations(ue,:)-AP_locations(ap,:)).^2)); 
        % if(UE_AP_dist(ap,ue)>D1)
        %     channelGaindB(ap,ue) = -1*constantTerm - 35*log10(UE_AP_dist(ap,ue));
        % elseif(UE_AP_dist(ap,ue)<=D1 && UE_AP_dist(ap,ue)>D0)
        %     channelGaindB(ap,ue) = -1*constantTerm - 15*log10(D1)-20*log10(UE_AP_dist(ap,ue));
        % else
        %     channelGaindB(ap,ue) = -1*constantTerm - 15*log10(D1)-20*log10(D0);
        % end
        channelGaindB(ap,ue) = -30.18-26*log10(UE_AP_dist(ap,ue));
        ricianFactor (ap,ue) = 10.^(1.3-0.003* UE_AP_dist(ap,ue));
    end
end
for ap1 = 1:M
    for ap2 = 1:M
        AP_AP_dist(ap1,ap2) = sqrt(sum((AP_locations(ap1,:)-AP_locations(ap2,:)).^2));
    end
end
for ue1 = 1:K
    for ue2 = 1:K
        UE_UE_dist(ue1,ue2) = sqrt(sum((UE_locations(ue1,:)-UE_locations(ue2,:)).^2));
    end
end

%decorr model
Cov_A = 2.^(-1*AP_AP_dist./d_DECORR);
Cov_B = 2.^(-1*UE_UE_dist./d_DECORR);

shadowing_a     = sqrtm(Cov_A)*(randn(M,1));
shadowing_b     = sqrtm(Cov_B)*(randn(K,1));
shadowing_z     = zeros(M,K);
for ue = 1:K
    for ap = 1:M
        %Generate shadow fading realizations
        shadowing_z(ap,ue)    = sigma_sf*(sqrt(delta)*shadowing_a(ap)+ sqrt(1-delta)*shadowing_b(ue));% correlated shadow fading model from paper
        % if (UE_AP_dist(ap,ue) > D0)
        channelGaindB(ap,ue) = channelGaindB(ap,ue) + shadowing_z(ap,ue);
        % end
    end
end
 channelGain_over_noise = 10.^(0.1*(channelGaindB-Noise_var_dB));
end