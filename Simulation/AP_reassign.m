function D_after_handoff =  AP_reassign(D, chgains, K_mmW, k_idx, p_fac, ue_idxs)
D_after_handoff = D;
M = size(chgains,1);
K = size(chgains,2);
ap_idxs = find(D(:,k_idx));
% for k = 1:(K-K_mmW)
loss_pc_arr = zeros(size(ue_idxs));
for kk_idx = 1:length(ue_idxs)
    k = ue_idxs(kk_idx);
    ap_idxs_k = find(D(:,k));
    ap_idxs_affected = intersect(ap_idxs,ap_idxs_k);
    loss_pc_arr(kk_idx) = 100*(sum(chgains(ap_idxs_affected,kk_idx))/sum(chgains(:,kk_idx)));
end
ue_idxs = ue_idxs(loss_pc_arr > 10);
for kk_idx = 1:length(ue_idxs)
    k = ue_idxs(kk_idx);
    ap_idxs_k = find(D(:,k));
    other_ap_idxs = setdiff(1:M,union(ap_idxs,ap_idxs_k));
    [~,other_ap_idxs_idxs] = sort(chgains(other_ap_idxs,k),'descend'); 
    other_ap_idxs = other_ap_idxs(other_ap_idxs_idxs);
    ap_idxs_affected = intersect(ap_idxs,ap_idxs_k);
    for m = 1:length(ap_idxs_affected)
        D_after_handoff(other_ap_idxs(m),k) = 1;
        if (p_fac == 1)
            D_after_handoff(ap_idxs_affected(m),k) = 0;
        end
    end
end
end