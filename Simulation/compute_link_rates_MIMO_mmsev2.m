function rate_dl = compute_link_rates_MIMO_mmsev2(p_d,D,M,Ntx,K_mmW,K,N_UE_mmW,N_UE_sub6,DeepMIMO_dataset_sub6,sub6ConnectionState,BW,NLoS_ue_idxs)
%% initialization of c
D_mmW_mmW = zeros(K_mmW,K_mmW,N_UE_mmW,N_UE_mmW);
D_mmW_sub6 = zeros(K_mmW,K-K_mmW,N_UE_mmW,N_UE_sub6);
D_sub6_mmW = zeros(K-K_mmW,K_mmW,N_UE_sub6,N_UE_mmW);
D_sub6_sub6 = zeros(K-K_mmW,K-K_mmW,N_UE_sub6,N_UE_sub6);
dl_mmse_precoder_mmW = zeros(M,K_mmW,Ntx,N_UE_mmW);
dl_mmse_precoder = zeros(M,K-K_mmW,Ntx,N_UE_sub6);
scaling_LP_mmse = zeros(M,K);
%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(BW) + noiseFigure;
noiseVariance = 10.^(noiseVariancedBm-30);
for m = 1:M
    for k = 1:K_mmW
        if (sub6ConnectionState(k) == 1)
            inv_matrix = noiseVariance*eye(Ntx);
            % inv_matrix = eye(Ntx);
            for q = 1:K_mmW
                if (D(m,q) == 1) && (sub6ConnectionState(q) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,"all")~=0)
                    % inv_matrix = inv_matrix + p_d*noiseVariance*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])';
                    inv_matrix = inv_matrix + p_d*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])';
                end
            end
            for q = 1:K-K_mmW
                if (D(m,q+K_mmW) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,"all")~=0)
                    % inv_matrix = inv_matrix + p_d*noiseVariance*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_sub6])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_sub6])';
                    inv_matrix = inv_matrix + p_d*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])';
                end
            end
            if (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,"all")~=0)
                dl_mmse_precoder_mmW(m,k,:,:) = reshape(dl_mmse_precoder_mmW(m,k,:,:),[Ntx,N_UE_mmW]) + p_d*inv_matrix\(reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,[Ntx,N_UE_mmW]));
            end
            if (D(m,k) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,"all")~=0)
                % scaling_LP_mmse(m,servedUEs) = scaling_LP_mmse(m,servedUEs) + norm(dl_mmse_precoder_mmW(m,k,:,:),'fro')^2;
                scaling_LP_mmse(m,k) = scaling_LP_mmse(m,k) + norm(dl_mmse_precoder_mmW(m,k,:,:),'fro')^2;
            end
        end
    end
    for k = 1:K-K_mmW
        inv_matrix = noiseVariance*eye(Ntx);
        % inv_matrix = eye(Ntx);
        for q = 1:K_mmW
            if (D(m,q) == 1) && (sub6ConnectionState(q) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,"all")~=0)
                % inv_matrix = inv_matrix + p_d*noiseVariance*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])';
                inv_matrix = inv_matrix + p_d*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])';
            end
        end
        for q = 1:K-K_mmW
            if (D(m,q+K_mmW) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,"all")~=0)
                % inv_matrix = inv_matrix +  p_d*noiseVariance*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_sub6])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_sub6])';
                inv_matrix = inv_matrix +  p_d*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])';
            end
        end
        if (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,"all")~=0)
            dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[Ntx,N_UE_sub6]) + p_d*inv_matrix\(reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,[Ntx,N_UE_sub6]));
        end
        if (D(m,k+K_mmW) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,"all")~=0)
            % scaling_LP_mmse(m,servedUEs) = scaling_LP_mmse(m,servedUEs) + norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
            scaling_LP_mmse(m,k+K_mmW) = scaling_LP_mmse(m,k+K_mmW) + norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
        end
    end
end
% dl_mmse_precoder_mmW = conj(channel_dl_mmW);
% dl_mmse_precoder = conj(channel_dl);
for m = 1:M
    for k = 1:K_mmW 
        if (D(m,k) == 1) && (sub6ConnectionState(k)==1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,"all")~=0)
            dl_mmse_precoder_mmW(m,k,:,:) = reshape(dl_mmse_precoder_mmW(m,k,:,:),[Ntx,N_UE_mmW])./sqrt(scaling_LP_mmse(m,k));
        end
    end
    for k = 1:K-K_mmW
        if (D(m,k+K_mmW) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,"all")~=0)
            dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[Ntx,N_UE_sub6])./sqrt(scaling_LP_mmse(m,k+K_mmW));
        end
    end
end
eta_eq = zeros(M,K);
if (K_mmW == 0) || all(sub6ConnectionState == 0)
    for m = 1:M
        term = 0;
        for k = 1+K_mmW:K
            if (D(m,k) == 1)
                % eta_eq(m,k) = 1/(Ntx*N_UE_sub6*sum(beta_uc(m,:)));
                term = term + trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:),[Ntx,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:),[Ntx,N_UE_sub6])');
            end
        end
        if (term > 0)
            eta_eq(m,:) = (1/term)*D(m,:);
        end
    end
else
    for m = 1:M
        term = 0;
        for k = 1:K
            if (D(m,k) == 1)
                if ((k<=K_mmW) && (sub6ConnectionState(k) == 1))
                    % eta_eq(m,k) = 1/(Ntx*(N_UE_mmW*(beta_uc(m,1:K_mmW).*sub6ConnectionState)+N_UE_sub6*sum(beta_uc(m,:))));
                    term = term + trace(reshape(dl_mmse_precoder_mmW(m,k,:,:),[Ntx,N_UE_mmW])*reshape(dl_mmse_precoder_mmW(m,k,:,:),[Ntx,N_UE_mmW])');
                elseif (k>K_mmW)
                    % eta_eq(m,k) = 1/(Ntx*(N_UE_mmW*(beta_uc(m,1:K_mmW).*sub6ConnectionState)+N_UE_sub6*sum(beta_uc(m,:))));
                    term = term + trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:),[Ntx,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:),[Ntx,N_UE_sub6])');
                end
            end
        end
        if (term > 0)
            eta_eq(m,:) = (1/term)*D(m,:);
            eta_eq(m,1:K_mmW) = eta_eq(m,1:K_mmW).*(sub6ConnectionState==1)';
        end
    end
end
for k = 1:K_mmW
    for q = 1:K_mmW
        for m = 1:M
            if (D(m,q) == 1)
                D_mmW_mmW(k,q,:,:) = reshape(D_mmW_mmW(k,q,:,:),[N_UE_mmW,N_UE_mmW]) + sqrt(eta_eq(m,q))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,[Ntx,N_UE_mmW])'*reshape(dl_mmse_precoder_mmW(m,q,:,:),[Ntx,N_UE_mmW]);
            end
        end
    end
    for q = 1:K-K_mmW
        for m = 1:M
            if (D(m,q+K_mmW) == 1)
                D_mmW_sub6(k,q,:,:) = reshape(D_mmW_sub6(k,q,:,:),[N_UE_mmW,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,[Ntx,N_UE_mmW])'*reshape(dl_mmse_precoder(m,q,:,:),[Ntx,N_UE_sub6]);
            end
        end
    end
end
for k = 1:K-K_mmW
    for q = 1:K_mmW
        for m = 1:M
            if (D(m,q) == 1)
                D_sub6_mmW(k,q,:,:) = reshape(D_sub6_mmW(k,q,:,:),[N_UE_sub6,N_UE_mmW]) + sqrt(eta_eq(m,q))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,[Ntx,N_UE_sub6])'*reshape(dl_mmse_precoder_mmW(m,q,:,:),[Ntx,N_UE_mmW]);
            end
        end
    end
    for q = 1:K-K_mmW
        for m = 1:M
            if (D(m,q+K_mmW) == 1)
                D_sub6_sub6(k,q,:,:) = reshape(D_sub6_sub6(k,q,:,:),[N_UE_sub6,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,[Ntx,N_UE_sub6])'*reshape(dl_mmse_precoder(m,q,:,:),[Ntx,N_UE_sub6]);
            end
        end
    end
end
DS_mmW = zeros(K_mmW,N_UE_mmW);
MSI_mmW = zeros(K_mmW,N_UE_mmW);
MUI_mmW = zeros(K_mmW,N_UE_mmW);
DS_sub6 = zeros(K-K_mmW,N_UE_sub6);
MSI_sub6 = zeros(K-K_mmW,N_UE_sub6);
MUI_sub6 = zeros(K-K_mmW,N_UE_sub6);

noise_mmW = abs(sqrt(0.5*noiseVariance)*(randn(K_mmW,N_UE_mmW) + 1j*randn(K_mmW,N_UE_mmW))).^2;
noise_sub6 = abs(sqrt(0.5*noiseVariance)*(randn(K-K_mmW,N_UE_sub6) + 1j*randn(K-K_mmW,N_UE_sub6))).^2;
snr_num_mmW = zeros(K_mmW,N_UE_mmW);
snr_den_mmW = zeros(K_mmW,N_UE_mmW);
snr_num_sub6 = zeros(K-K_mmW,N_UE_sub6);
snr_den_sub6 = zeros(K-K_mmW,N_UE_sub6);
rate_dl = zeros(K,1);
for k = 1:K_mmW
%     if (sub6ConnectionState(k)==1 || k==ue_idx)
    if (sub6ConnectionState(k)==1)
        for n = 1:N_UE_mmW
%             DS_mmW(k,n) = p_d*norm(reshape(D_mmW_mmW(k,k,n,:),[1,N_UE_mmW]))^2;
            DS_mmW(k,n) = p_d*(abs(D_mmW_mmW(k,k,n,n)))^2;
            for nn = 1:N_UE_mmW
%                 if (nn~=n)
                if (abs(D_mmW_mmW(k,k,nn,nn))<abs(D_mmW_mmW(k,k,n,n)))
%                     MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,k,nn,:),[1,N_UE_mmW]))^2;
                    MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*(abs(D_mmW_mmW(k,k,n,nn)))^2;
                end
            end
            for q = 1:K_mmW
                if (q~=k && sub6ConnectionState(q)==1)
                  MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,q,n,:),[1,N_UE_mmW]))^2;
                end
            end
            for q = 1:K-K_mmW
               MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_sub6(k,q,n,:),[1,N_UE_sub6]))^2;
            end
            snr_num_mmW(k,n) = DS_mmW(k,n);
            snr_den_mmW(k,n) = MSI_mmW(k,n) + MUI_mmW(k,n) + noise_mmW(k,n);
            rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
        end
    end
end
for k = 1:K-K_mmW
    for n = 1:N_UE_sub6
%         DS_sub6(k,n) = p_d*norm(reshape(D_sub6_sub6(k,k,n,:),[1,N_UE_sub6]))^2;
        DS_sub6(k,n) = p_d*(abs(D_sub6_sub6(k,k,n,n)))^2;
        for nn = 1:N_UE_sub6
%             if (nn~=n)
%             if (nn>n)
            if (abs(D_sub6_sub6(k,k,nn,nn))<abs(D_sub6_sub6(k,k,n,n)))
%                 MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,k,nn,:),[1,N_UE_sub6]))^2;
                MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*(abs(D_sub6_sub6(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_mmW
            if (q~=k && sub6ConnectionState(q)==1)
              MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_mmW(k,q,n,:),[1,N_UE_mmW]))^2;
            end
        end
        for q = 1:K-K_mmW
           MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,q,n,:),[1,N_UE_sub6]))^2;
        end
        snr_num_sub6(k,n) = DS_sub6(k,n);
        snr_den_sub6(k,n) = MSI_sub6(k,n) + MUI_sub6(k,n) + noise_sub6(k,n);
        rate_dl(k+K_mmW) = rate_dl(k+K_mmW) + BW*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
    end
end
end