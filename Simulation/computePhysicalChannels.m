%function [BS_H_fr_dl,  BS_H_fr_ul]= computePhysicalChannels(params)
function [BS_H_fr_dl_mmw, BS_H_fr_dl_sub6]= computePhysicalChannels(params)

UE_locations = params.UE_locations;
ue_height = params.hr;
BS_height = params.ht;
locations_BS = params.locationsBS;
CenterFreq = 30e9;
BS_Nt = 256;
M = params.numGNB;
ue_Nt = 32;
ue_Nr = 32;
K = params.num_sc_mmw; % no of subcarriers
% scs = 120e3; %sub carrier spacing
scs = params.scs_mmw;
BW = K*scs;
Kd = params.numUE;
if (Kd > 0)
    [~, ~, BS_H_fr_dl_mmw, ~] = sample_generate_3gpp_channel('3GPP_38.901_UMa_LOS', CenterFreq, BW, M, Kd, BS_height, ue_height, BS_Nt, ue_Nr, K, UE_locations, locations_BS);
end
BS_height = params.ht_sub6;
locations_BS = params.locationsBS_sub6;
CenterFreq = 2.5e9;
BS_Nt = 256;
M = params.numGNB_sub6;
ue_Nt = 32;
ue_Nr = 32;
K = params.num_sc_sub6; % no of subcarriers
% scs = 15e3; %sub carrier spacing
scs = params.scs_sub6;
BW = K*scs;
Kd = params.numUE;
if (Kd > 0)
    [~, ~, BS_H_fr_dl_sub6, ~] = sample_generate_3gpp_channel('3GPP_38.901_UMa_NLOS', CenterFreq, BW, M, Kd, BS_height, ue_height, BS_Nt, ue_Nr, K, UE_locations, locations_BS);
end
% % BS_Nr = 256;
% % ue_Nt = 32;
% K = 792; % no of subcarriers
% scs = 120e3; %sub carrier spacing
% BW = K*scs;
% max_BS_ue_dist = 100;
% min_BS_ue_dist = 10;
% if (Ku > 0)
%     [~, ~, BS_H_fr_ul, c_obj_ul] = sample_generate_3gpp_channel('3GPP_38.901_UMa_LOS', CenterFreq, BW, Ku, BS_height, ue_height, Nrx, ue_Nt, K, max_BS_ue_dist, min_BS_ue_dist);
% end