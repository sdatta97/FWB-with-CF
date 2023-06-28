function link_rates = compute_link_rates(ch_coeffs, scs, idx_update, ue_idx)
M = size(ch_coeffs,2);
link_rates = zeros(M,1);
% scs = 120e3; %sub carrier spacing
% N0 = 1;
N0_dBm = -228.6 + 30 + 10*log10(290);
kappa_dB = 9;
kappa = 10^(0.1*kappa_dB);
N0 = 10^(0.1*N0_dBm)*10^(-3);
Pt_dBm = -10;
Pt = 10^(0.1*Pt_dBm)*10^(-3);
Nr = size(ch_coeffs{ue_idx,1},1);
Nt = size(ch_coeffs{ue_idx,1},2);
Nc = size(ch_coeffs{ue_idx,1},3);
BW = Nc*scs;
for i = 1:M
    for j = 1:Nc 
        chgains = svd(ch_coeffs{ue_idx,i}(:,:,j));
        % if idx_update>0
        if idx_update~=i
            % chgains = sqrt(1/(1+kappa)).*chgains;
            link_rates(i) = link_rates(i) + scs*log2(1+Pt*sum(chgains.^2)/N0*scs);
        end
    end
end
end