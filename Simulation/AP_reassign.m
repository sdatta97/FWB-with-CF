function D_after_handoff =  AP_reassign(D, chgains, K_mmW, k_idx, p_fac)
D_after_handoff = D;
M = size(chgains,1);
K = size(chgains,2);
ap_idxs = find(D(:,k_idx));
for k = 1:(K-K_mmW)
    ap_idxs_k = find(D(:,k+K_mmW));
    other_ap_idxs = setdiff(1:M,union(ap_idxs,ap_idxs_k));
    [~,other_ap_idxs_idxs] = sort(chgains(other_ap_idxs,k+K_mmW),'descend'); 
    other_ap_idxs = other_ap_idxs(other_ap_idxs_idxs);
    ap_idxs_affected = intersect(ap_idxs,ap_idxs_k);
    for m = 1:length(ap_idxs_affected)
        D_after_handoff(other_ap_idxs(m),k+K_mmW) = 1;
        if (p_fac == 1)
            D_after_handoff(ap_idxs_affected(m),k+K_mmW) = 0;
        end
    end
end
end