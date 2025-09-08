clear
% load("../../Desktop/data-mmW.mat","DeepMIMO_dataset")
% DeepMIMO_dataset_mmW = DeepMIMO_dataset;
% bs_dataset = DeepMIMO_dataset{2};
% bs_loc = bs_dataset.loc;
% user_locs_mmW = zeros(length(bs_dataset.user), 2);
% user_los_mmW = zeros(length(bs_dataset.user), 1);
% for j = 1:length(bs_dataset.user)
%     user_locs_mmW(j, :) = bs_dataset.user{j}.loc(1:2);
%     user_los_mmW(j) = bs_dataset.user{j}.LoS_status;
% end
% clear DeepMIMO_dataset
load("../../Desktop/data-sub6.mat","DeepMIMO_dataset")
DeepMIMO_dataset_sub6 = DeepMIMO_dataset;
bs_dataset = DeepMIMO_dataset{3};
user_locs_sub6 = zeros(length(bs_dataset.user), 2);
user_los_sub6 = zeros(length(bs_dataset.user), 1);
for j = 1:length(bs_dataset.user)
    user_locs_sub6(j, :) = bs_dataset.user{j}.loc(1:2);
    user_los_sub6(j) = bs_dataset.user{j}.LoS_status;
end

clear DeepMIMO_dataset
p_fac = 10;
loss_arr = 0:0.1:1;
M = 12;
K = 50;
K_mmW = 5;
r_min_sub6 = 35e6;
BW = 1e8;
TAU_FAC = 1;
Ntx = 64;
N_UE_mmW = 1;
N_UE_sub6 = 1;
p_d = 10^(3.6)*Ntx; % params.rho_tot; % 1*K;
D = ones(M,K); %randi([0,1],[M,K]); % params.D
Lmax = 4;
all_BS_loc = zeros(M,2);
for m = 1:M
    all_BS_loc(m,:) = DeepMIMO_dataset_sub6{m}.loc(1:2);
end
num_samples = 10000;
rate_dl_mmW = zeros([num_samples,size(loss_arr,2)]);
rate_dl_sub6 = zeros([num_samples, size(loss_arr,2)]);
for i = 1:num_samples
    NLoS_ue_idxs = find(user_los_sub6==0);
    NLoS_ue_idxs = NLoS_ue_idxs(randsample(numel(NLoS_ue_idxs),K));
    user_locs_sub6 = user_locs_sub6(NLoS_ue_idxs,:);
    user_los_sub6 = user_los_sub6(NLoS_ue_idxs);
    
    % [dist_loc_sub6, NLoS_ue_idxs] = mink(sqrt(sum((user_locs_sub6-repmat(bs_dataset.loc(1:2),[size(user_locs_sub6,1),1])).^2,2)),K);
    dist_loc_sub6 = sqrt(sum((user_locs_sub6-repmat(bs_dataset.loc(1:2),[size(user_locs_sub6,1),1])).^2,2));
    for k = 1:K
        [~, idxs] = sort(sum((all_BS_loc - repmat(user_locs_sub6(k,:),[M,1])).^2,2));
        idxs_not_chosen = idxs((Lmax+1):end);
        D(idxs_not_chosen,k) = 0;
    end
    
    ue_idx = 1;
    sub6ConnectionState = zeros(K_mmW,1);
    sub6ConnectionState(ue_idx) = 1;
    ap_idxs = find(D(:,ue_idx));
    ue_idxs = [];
    for j = 1:length(ap_idxs)
        ap = ap_idxs(j);
        ue_idxs = union(ue_idxs, find(D(ap,:)));
    end
    ue_rearranged = ue_idxs ((1+K_mmW):end);
    for j=1:length(loss_arr)
        loss_pc_thresh = loss_arr(j);
        rate_dl_before_handoff = compute_link_rates_MIMO_mmsev2(p_d,D,M,Ntx,K_mmW,K,N_UE_mmW,N_UE_sub6,DeepMIMO_dataset_sub6,zeros(K_mmW,1),BW,NLoS_ue_idxs);
        lb = quantile(rate_dl_before_handoff(ue_rearranged)./BW,loss_pc_thresh);
        bw_alloc = BW - r_min_sub6/lb;
        if (bw_alloc < 0) %|| isnan(bw_alloc)
            bw_alloc = 0;
            p_fac = 1;
            rate_dl_after_handoff = rate_dl_before_handoff;
        elseif isnan(bw_alloc)
            rate_dl_after_handoff = compute_link_rates_MIMO_mmsev2(p_d,D,M,Ntx,K_mmW,K,N_UE_mmW,N_UE_sub6,DeepMIMO_dataset_sub6,sub6ConnectionState,BW,NLoS_ue_idxs); 
            lb = quantile(rate_dl_after_handoff((1+K_mmW):end),loss_pc_thresh);                                
        else 
            ues_not_affected = setdiff((1+K_mmW):K,ue_rearranged);
            num_sc_sub6 = 2;
            scs_sub6 = zeros(num_sc_sub6,1);
            scs_sub6(1) = bw_alloc;
            scs_sub6(2) = BW - bw_alloc;
            user_sc_alloc = zeros(K,1);     
            user_sc_alloc(find(sub6ConnectionState),1) = 1;
            user_sc_alloc(find(sub6ConnectionState),2) = 0;
            user_sc_alloc(ues_not_affected,1) = 1;
            user_sc_alloc(ues_not_affected,2) = 1;
            user_sc_alloc(ue_rearranged,1) = 0;
            user_sc_alloc(ue_rearranged,2) = 1;
            ues_sharing = union(((1:K_mmW).*sub6ConnectionState),ues_not_affected);
            rate_dl_after_handoff = compute_link_rates_MIMOv4v2(p_d,p_fac,D,M,Ntx,K_mmW,K,N_UE_mmW,N_UE_sub6,DeepMIMO_dataset_sub6,sub6ConnectionState,user_sc_alloc,scs_sub6,num_sc_sub6,NLoS_ue_idxs);    
        end
        rate_dl_mmW(i,j) = rate_dl_after_handoff(ue_idx);
        rate_dl_sub6(i,j) = median(rate_dl_after_handoff(1+K_mmW:end));
    end
end