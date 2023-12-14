clear all;

CenterFreq = 30e9;
num_ue = 1;
BS_height = 10; % in meters
ue_height = 1.5;
BS_Nt = 256;
ue_Nr = 32;
K = 792; % no of subcarriers
scs = 120e3; %sub carrier spacing
BW = K*scs;
max_BS_ue_dist = 100;
min_BS_ue_dist = 10;
[~, ~, BS_H_fr, c_obj] = sample_generate_3gpp_channel('3GPP_38.901_UMi_LOS', CenterFreq, BW, num_ue, BS_height, ue_height, BS_Nt, ue_Nr, K, max_BS_ue_dist, min_BS_ue_dist);