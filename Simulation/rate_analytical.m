function rate_dl = rate_analytical(params)
M = params.numGNB_sub6;
K_mmW = params.numUE;
K = params.numUE + params.numUE_sub6;
TAU_FAC = params.preLogFactor;
N_AP = params.num_antennas_per_gNB;
N_UE_mmW = params.N_UE_mmW;
N_UE_sub6 = params.N_UE_sub6;
p_d = params.rho_tot; % 1*K;
p_fac = params.p_fac;
p_fac_rearrange = params.p_fac_rearrange;
D = params.D;
R_gNB = params.R_gNB;

%Prepare cell to store the AP indices serving a specific UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specific UE
NoServ = cell(K,1);
%Construc the above array and cells
for k = 1:K
    servingAPs = find(D(:,k)==1);
    NoservingAPs = find(D(:,k)==0);
    
    Serv{k} = servingAPs;
    NoServ{k} = NoservingAPs;   
end

rate_dl = zeros(K,1);
DS_mmW = zeros(K_mmW,N_UE_mmW);
MSI_mmW = zeros(K_mmW,N_UE_mmW);
MUI_mmW = zeros(K_mmW,N_UE_mmW);
DS_sub6 = zeros(K-K_mmW,N_UE_sub6);
MSI_sub6 = zeros(K-K_mmW,N_UE_sub6);
MUI_sub6 = zeros(K-K_mmW,N_UE_sub6);   
noise_mmW = abs(sqrt(0.5)*(randn(K_mmW,N_UE_mmW) + 1j*randn(K_mmW,N_UE_mmW))).^2;
noise_sub6 = abs(sqrt(0.5)*(randn(K-K_mmW,N_UE_sub6) + 1j*randn(K-K_mmW,N_UE_sub6))).^2;
snr_num_mmW = zeros(K_mmW,N_UE_mmW);
snr_den_mmW = zeros(K_mmW,N_UE_mmW);
snr_num_sub6 = zeros(K-K_mmW,N_UE_sub6);
snr_den_sub6 = zeros(K-K_mmW,N_UE_sub6);
for k = 1:K_mmW
    if (sub6ConnectionState(k)==1) 
        for n = 1:N_UE_mmW              
            DS_mmW(k,n) = p_d*(abs(D_mmW_mmW(k,k,n,n)))^2;
            for nn = 1:N_UE_mmW
               MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*(abs(D_mmW_mmW(k,k,n,nn)))^2;
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
            rate_dl(k) = rate_dl(k) + TAU_FAC*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
        end
    end
end
for k = 1:K-K_mmW
    for n = 1:N_UE_sub6
        DS_sub6(k,n) = p_d*(abs(D_sub6_sub6(k,k,n,n)))^2;
        for nn = 1:N_UE_sub6
           MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*(abs(D_sub6_sub6(k,k,n,nn)))^2;
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
        rate_dl(k+K_mmW) = rate_dl(k+K_mmW) + TAU_FAC*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
    end
end
end