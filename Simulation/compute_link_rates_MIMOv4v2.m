function rate_dl = compute_link_rates_MIMOv4v2(p_d,p_fac,D,M,Ntx,K_mmW,K,N_UE_mmW,N_UE_sub6,DeepMIMO_dataset_sub6,sub6ConnectionState,user_sc_alloc,scs,num_sc_sub6,NLoS_ue_idxs)
dl_mmse_precoder_mmW = zeros([M,K_mmW,Ntx,N_UE_mmW,num_sc_sub6]);
dl_mmse_precoder = zeros([M,K-K_mmW,Ntx,N_UE_sub6,num_sc_sub6]);
scaling_LP_mmse = zeros(M,K,num_sc_sub6);
%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(scs) + noiseFigure;
noiseVariance = 10.^(noiseVariancedBm-30);
for m = 1:M
    for n = 1:num_sc_sub6
        for k = 1:K_mmW
            if (sub6ConnectionState(k) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,"all")~=0) && (user_sc_alloc(k,n) == 1)
                inv_matrix = noiseVariance(n)*eye(Ntx);
                %inv_matrix = eye(Ntx);
                for q = 1:K_mmW
                    if (D(m,q)==1) && (sub6ConnectionState(q) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,"all")~=0) && (user_sc_alloc(q,n) == 1)
                        % inv_matrix = inv_matrix + p_d*noiseVariance(n)*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])';
                        inv_matrix = inv_matrix + p_d*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])';
                    end
                end
                for q = 1:K-K_mmW
                    if (D(m,q+K_mmW) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,"all")~=0) && (user_sc_alloc(q+K_mmW,n) == 1)
                        % inv_matrix = inv_matrix + p_d*noiseVariance(n)*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])';
                        inv_matrix = inv_matrix + p_d*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])';
                    end
                end
                if (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,"all")~=0)
                    dl_mmse_precoder_mmW(m,k,:,:,n) = reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[Ntx,N_UE_mmW]) + p_d*inv_matrix\(reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,[Ntx,N_UE_mmW]));
                end
                if (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,"all")~=0)
                    % scaling_LP_mmse(m,servedUEs) = scaling_LP_mmse(m,servedUEs) + norm(dl_mmse_precoder_mmW(m,k,:,:),'fro')^2;
                    scaling_LP_mmse(m,k,n) = scaling_LP_mmse(m,k,n) + norm(dl_mmse_precoder_mmW(m,k,:,:,n),'fro')^2;
                end
            end
        end
        for k = 1:K-K_mmW
            if (user_sc_alloc(k+K_mmW,n) == 1) 
                inv_matrix = noiseVariance(n)*eye(Ntx);
                % inv_matrix = eye(Ntx);
                for q = 1:K_mmW
                    if (D(m,q) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,"all")~=0) && (sub6ConnectionState(q) == 1) && (user_sc_alloc(q,n) == 1)
                        % inv_matrix = inv_matrix + p_d*noiseVariance(n)*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])';
                        inv_matrix = inv_matrix + p_d*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q)}.channel,[Ntx,N_UE_mmW])';
                    end
                end
                for q = 1:K-K_mmW
                    if (D(m,q+K_mmW) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,"all")~=0) && (user_sc_alloc(q+K_mmW,n) == 1)
                        % inv_matrix = inv_matrix +  p_d*noiseVariance(n)*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])';
                        inv_matrix = inv_matrix +  p_d*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(q+K_mmW)}.channel,[Ntx,N_UE_sub6])';
                    end
                end
                if (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,"all")~=0)
                    dl_mmse_precoder(m,k,:,:,n) = reshape(dl_mmse_precoder(m,k,:,:,n),[Ntx,N_UE_sub6]) + p_d*inv_matrix\(reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,[Ntx,N_UE_sub6]));
                end
                if (D(m,k+K_mmW) == 1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,"all")~=0)
                    % scaling_LP_mmse(m,servedUEs) = scaling_LP_mmse(m,servedUEs) + norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
                    scaling_LP_mmse(m,k+K_mmW,n) = scaling_LP_mmse(m,k+K_mmW,n) + norm(dl_mmse_precoder(m,k,:,:,n),'fro')^2;
                end
            end
        end
    end
end
% dl_mmse_precoder_mmW = conj(channel_est_dl_mmW);
% dl_mmse_precoder = conj(channel_est_dl);
for m = 1:M
    for n = 1:num_sc_sub6
        for k = 1:K_mmW
            if (D(m,k)==1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,"all")~=0) && (user_sc_alloc(k,n) == 1) && (sub6ConnectionState(k)==1)
                dl_mmse_precoder_mmW(m,k,:,:,n) = reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[Ntx,N_UE_mmW])./sqrt(scaling_LP_mmse(m,k,n));
            end
        end
        for k = 1:K-K_mmW
            if (D(m,k+K_mmW)==1) && (sum(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,"all")~=0) && (user_sc_alloc(k+K_mmW,n) == 1)
                dl_mmse_precoder(m,k,:,:,n) = reshape(dl_mmse_precoder(m,k,:,:,n),[Ntx,N_UE_sub6])./sqrt(scaling_LP_mmse(m,k+K_mmW,n));
            end
        end
    end
end

%% initialization of c
eta_eq = zeros(M,K,num_sc_sub6);
if (K_mmW == 0) || all(sub6ConnectionState == 0)
    for m = 1:M
        for n = 1:num_sc_sub6
            term = 0;
            for k = 1+K_mmW:K
                if (D(m,k)==1) && (user_sc_alloc(k,n) == 1)
                        %term = (Ntx*N_UE_sub6*num_sc_sub6*beta_uc(m,:)*user_sc_alloc(:,n));
                        % term = num_sc_sub6*Ntx*(N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                      
                        % term = term + trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[Ntx,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[Ntx,N_UE_sub6])');
                        term = term + num_sc_sub6*trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[Ntx,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[Ntx,N_UE_sub6])');
                end
            end
            if (term > 0)
                eta_eq(m,:,n) = (1/term)*(D(m,:)'.*(user_sc_alloc(:,n)==1));
            end
        end
    end
else
    % ues_not_rearranged = setdiff((1+K_mmW):K,ue_rearranged);
    for m = 1:M
        for n = 1:num_sc_sub6
            term = 0;
            for k = 1:K
                if (D(m,k)==1) && (user_sc_alloc(k,n) == 1)
                    if (k<=K_mmW) && (sub6ConnectionState(k) == 1) 
    %                         term = (Ntx*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*beta_uc(m,(1+K_mmW):K)*user_sc_alloc((1+K_mmW):K,n)));
                        % term = num_sc_sub6*Ntx*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                                              
                        % term = term + p_fac*trace(reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[Ntx,N_UE_mmW])*reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[Ntx,N_UE_mmW])');
                        term = term + num_sc_sub6*p_fac*trace(reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[Ntx,N_UE_mmW])*reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[Ntx,N_UE_mmW])');
                    elseif (k>K_mmW)
%                     eta_eq(m,k) = 1./(Ntx*(N_UE_mmW*p_fac*beta_uc(m,1:K_mmW)+N_UE_sub6*sum(beta_uc(m,2:K))));
%                         term = (Ntx*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*beta_uc(m,(1+K_mmW):K)*user_sc_alloc((1+K_mmW):K,n)));
                        % term = num_sc_sub6*Ntx*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                        
                        % term = term + trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[Ntx,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[Ntx,N_UE_sub6])');
                        term = term + num_sc_sub6*trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[Ntx,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[Ntx,N_UE_sub6])');
                    end
                end
            end
            if (term > 0)
                eta_eq(m,:,n) = (1/term)*(D(m,:)'.*(user_sc_alloc(:,n)==1));
                eta_eq(m,1:K_mmW,n) = p_fac*(eta_eq(m,1:K_mmW,n).*(sub6ConnectionState==1)');
            end
        end
    end
end
D_mmW_mmW = zeros(K_mmW,K_mmW,N_UE_mmW,N_UE_mmW,num_sc_sub6);
D_mmW_sub6 = zeros(K_mmW,K-K_mmW,N_UE_mmW,N_UE_sub6,num_sc_sub6);
D_sub6_mmW = zeros(K-K_mmW,K_mmW,N_UE_sub6,N_UE_mmW,num_sc_sub6);
D_sub6_sub6 = zeros(K-K_mmW,K-K_mmW,N_UE_sub6,N_UE_sub6,num_sc_sub6);
for n_idx = 1:num_sc_sub6
    for k = 1:K_mmW
        for q = 1:K_mmW
            for m = 1:M
                if (D(m,q)==1) && (sub6ConnectionState(q) == 1) && (user_sc_alloc(q,n_idx) == 1)
                    % D_mmW_mmW(k,q,:,:,n_idx) = reshape(D_mmW_mmW(k,q,:,:,n_idx),[N_UE_mmW,N_UE_mmW]) + sqrt(eta_eq(m,q,n_idx))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,[Ntx,N_UE_mmW]).'*reshape(conj(channel_est_dl_mmW(m,q,:,:)),[Ntx,N_UE_mmW]);
                    D_mmW_mmW(k,q,:,:,n_idx) = reshape(D_mmW_mmW(k,q,:,:,n_idx),[N_UE_mmW,N_UE_mmW]) + sqrt(eta_eq(m,q,n_idx))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,[Ntx,N_UE_mmW])'*reshape(dl_mmse_precoder_mmW(m,q,:,:,n_idx),[Ntx,N_UE_mmW]);
                end
            end
        end
        for q = 1:K-K_mmW
            for m = 1:M
                if (D(m,q+K_mmW)==1) && (user_sc_alloc(q+K_mmW,n_idx) == 1)
                    % D_mmW_sub6(k,q,:,:,n_idx) = reshape(D_mmW_sub6(k,q,:,:,n_idx),[N_UE_mmW,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW,n_idx))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,[Ntx,N_UE_mmW]).'*reshape(conj(channel_est_dl(m,q,:,:)),[Ntx,N_UE_sub6]);
                    D_mmW_sub6(k,q,:,:,n_idx) = reshape(D_mmW_sub6(k,q,:,:,n_idx),[N_UE_mmW,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW,n_idx))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k)}.channel,[Ntx,N_UE_mmW])'*reshape(dl_mmse_precoder(m,q,:,:,n_idx),[Ntx,N_UE_sub6]);
                end
            end
        end
    end
    for k = 1:K-K_mmW
        for q = 1:K_mmW
            for m = 1:M
                if (D(m,q)==1) && (sub6ConnectionState(q) == 1) && (user_sc_alloc(q,n_idx) == 1)
                    % D_sub6_mmW(k,q,:,:,n_idx) = reshape(D_sub6_mmW(k,q,:,:,n_idx),[N_UE_sub6,N_UE_mmW]) + sqrt(eta_eq(m,q,n_idx))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,[Ntx,N_UE_sub6]).'*reshape(conj(channel_est_dl_mmW(m,q,:,:)),[Ntx,N_UE_mmW]);
                    D_sub6_mmW(k,q,:,:,n_idx) = reshape(D_sub6_mmW(k,q,:,:,n_idx),[N_UE_sub6,N_UE_mmW]) + sqrt(eta_eq(m,q,n_idx))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,[Ntx,N_UE_sub6])'*reshape(dl_mmse_precoder_mmW(m,q,:,:,n_idx),[Ntx,N_UE_mmW]);
                end
            end
        end
        for q = 1:K-K_mmW
            for m = 1:M
                if (D(m,q+K_mmW)==1) && (user_sc_alloc(q+K_mmW,n_idx) == 1)
                    % D_sub6_sub6(k,q,:,:,n_idx) = reshape(D_sub6_sub6(k,q,:,:,n_idx),[N_UE_sub6,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW,n_idx))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,[Ntx,N_UE_sub6]).'*reshape(conj(channel_est_dl(m,q,:,:)),[Ntx,N_UE_sub6]);
                    D_sub6_sub6(k,q,:,:,n_idx) = reshape(D_sub6_sub6(k,q,:,:,n_idx),[N_UE_sub6,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW,n_idx))*reshape(DeepMIMO_dataset_sub6{m}.user{NLoS_ue_idxs(k+K_mmW)}.channel,[Ntx,N_UE_sub6])'*reshape(dl_mmse_precoder(m,q,:,:,n_idx),[Ntx,N_UE_sub6]);
                end
            end
        end
    end
end
rate_dl = zeros(K,1);
for n_idx = 1:num_sc_sub6
    DS_mmW = zeros(K_mmW,N_UE_mmW);
    MSI_mmW = zeros(K_mmW,N_UE_mmW);
    MUI_mmW = zeros(K_mmW,N_UE_mmW);
    DS_sub6 = zeros(K-K_mmW,N_UE_sub6);
    MSI_sub6 = zeros(K-K_mmW,N_UE_sub6);
    MUI_sub6 = zeros(K-K_mmW,N_UE_sub6);
    
    noise_mmW = abs(sqrt(0.5*noiseVariance(n_idx))*(randn(K_mmW,N_UE_mmW) + 1j*randn(K_mmW,N_UE_mmW))).^2;
    noise_sub6 = abs(sqrt(0.5*noiseVariance(n_idx))*(randn(K-K_mmW,N_UE_sub6) + 1j*randn(K-K_mmW,N_UE_sub6))).^2;
    snr_num_mmW = zeros(K_mmW,N_UE_mmW);
    snr_den_mmW = zeros(K_mmW,N_UE_mmW);
    snr_num_sub6 = zeros(K-K_mmW,N_UE_sub6);
    snr_den_sub6 = zeros(K-K_mmW,N_UE_sub6);
    for k = 1:K_mmW
    %     if (sub6ConnectionState(k)==1 || k==ue_idx)
        if (sub6ConnectionState(k)==1) && (user_sc_alloc(k,n_idx) == 1)
            for n = 1:N_UE_mmW              
    %             DS_mmW(k,n) = p_d*norm(reshape(D_mmW_mmW(k,k,n,:),[1,N_UE_mmW]))^2;
                DS_mmW(k,n) = p_d*(abs(D_mmW_mmW(k,k,n,n,n_idx)))^2;
                for nn = 1:N_UE_mmW
    %                 if (nn~=n)
                    if (abs(D_mmW_mmW(k,k,nn,nn,n_idx))<abs(D_mmW_mmW(k,k,n,n,n_idx)))
    %                     MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,k,nn,:),[1,N_UE_mmW]))^2;
                        MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*(abs(D_mmW_mmW(k,k,n,nn,n_idx)))^2;
                    end
                end
                for q = 1:K_mmW
                    if (q~=k && sub6ConnectionState(q)==1)
                      MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,q,n,:,n_idx),[1,N_UE_mmW]))^2;
                    end
                end
                for q = 1:K-K_mmW
                   MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_sub6(k,q,n,:,n_idx),[1,N_UE_sub6]))^2;
                end
                snr_num_mmW(k,n) = DS_mmW(k,n);
                snr_den_mmW(k,n) = MSI_mmW(k,n) + MUI_mmW(k,n) + noise_mmW(k,n);
                rate_dl(k) = rate_dl(k) + scs(n_idx)*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
            end
        end
    end
    for k = 1:K-K_mmW
        for n = 1:N_UE_sub6
           if (user_sc_alloc(k+K_mmW,n_idx) == 1)
    
        %         DS_sub6(k,n) = p_d*norm(reshape(D_sub6_sub6(k,k,n,:),[1,N_UE_sub6]))^2;
                DS_sub6(k,n) = p_d*(abs(D_sub6_sub6(k,k,n,n,n_idx)))^2;
                for nn = 1:N_UE_sub6
        %             if (nn~=n)
        %             if (nn>n)
                    if (abs(D_sub6_sub6(k,k,nn,nn,n_idx))<abs(D_sub6_sub6(k,k,n,n,n_idx)))
        %                 MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,k,nn,:),[1,N_UE_sub6]))^2;
                        MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*(abs(D_sub6_sub6(k,k,n,nn,n_idx)))^2;
                    end
                end
                for q = 1:K_mmW
                    if (q~=k && sub6ConnectionState(q)==1)
                      MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_mmW(k,q,n,:,n_idx),[1,N_UE_mmW]))^2;
                    end
                end
                for q = 1:K-K_mmW
                   MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,q,n,:,n_idx),[1,N_UE_sub6]))^2;
                end
                snr_num_sub6(k,n) = DS_sub6(k,n);
                snr_den_sub6(k,n) = MSI_sub6(k,n) + MUI_sub6(k,n) + noise_sub6(k,n);
                rate_dl(k+K_mmW) = rate_dl(k+K_mmW) + scs(n_idx)*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
            end
        end
    end
end
end