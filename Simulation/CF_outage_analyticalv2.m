function p_out = CF_outage_analyticalv2(params,ue_idx,lambda_BS,lambda_UE)
    p_out = 0.5;
    D = params.D;
    ue_rearranged = params.ue_rearranged;
    p_fac = params.p_fac;
    p_fac_rearrange = params.p_fac_rearrange;
    num_sc_sub6 = params.num_sc_sub6;
    scs_sub6 = params.scs_sub6;
    N = params.num_antennas_per_gNB;
    N_UE_sub6 = params.N_UE_sub6;
    N_UE_mmW = params.N_UE_mmW;
    M = params.numGNB_sub6;
    K = params.numUE+params.numUE_sub6;
    K_mmW = params.numUE;
    p_d = params.rho_tot;
    B = params.Band;
    rmin = params.r_min_sub6 (1);
    BETA = params.BETA;
    beta_uc = zeros(size(BETA));
    sub6ConnectionState = params.sub6ConnectionState;
    %Prepare array to store the number of APs serving a specficic UE
    La = zeros(K,1);
    %Prepare cell to store the AP indices serving a specficic UE
    Serv = cell(K,1);
    %Prepare cell to store the AP indices not serving a specficic UE
    NoServ = cell(K,1);
    %Construc the above array and cells
    for k = 1:K
        servingAPs = find(D(:,k)==1);
        NoservingAPs = find(D(:,k)==0);
        
        Serv{k} = servingAPs;
        NoServ{k} = NoservingAPs;
        
        La(k) = length(servingAPs);
        beta_uc(:,k) = BETA(:,k).*D(:,k);
    end
    
    %% initialization of c
    eta_eq = zeros(M,K,num_sc_sub6);
    user_sc_alloc = params.user_sc_alloc;
    N_AP = params.num_antennas_per_gNB;
    N_tilde_num = zeros(K,1);
    N_tilde_den = zeros(K,1);
    beta_tilde_num = zeros(K,1);
    beta_tilde_den = zeros(K,1);
    k_tilde_num = zeros(K,1);
    k_tilde_den = zeros(K,1);
    theta_tilde_num = zeros(K,1);
    theta_tilde_den = zeros(K,1);    
    if ((K_mmW == 0) || (sub6ConnectionState == zeros(K_mmW,1)))
        ues_not_rearranged = setdiff((1+K_mmW):K,ue_rearranged);
        for m = 1:M
            for k = 1+K_mmW:K
                if ismember(m,Serv{k})
                    for n = 1:num_sc_sub6
                        %term = (N_AP*N_UE_sub6*num_sc_sub6*beta_uc(m,:)*user_sc_alloc(:,n));
                        term = num_sc_sub6*N_AP*(N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                      
                        if term > 0
                            if ismember(k,ue_rearranged) 
                                eta_eq(m,k,n) = p_fac_rearrange/term;
                            else
                                eta_eq(m,k,n) = 1/term;
                            end
                        end
                        N_tilde_num(k) = N_tilde_num(k) + sqrt(p_d*eta_eq(m,k,n))*beta_uc(m,k);
                        N_tilde_den(k) = N_tilde_den(k) + p_d*eta_eq(m,k,n)*(beta_uc(m,k))^2;
                        beta_tilde_num(k) = beta_tilde_num(k) + p_d*eta_eq(m,k,n)*(beta_uc(m,k))^2;
                        beta_tilde_den(k) = beta_tilde_den(k) + sqrt(p_d*eta_eq(m,k,n))*beta_uc(m,k);
                        k_tilde_num(k) = k_tilde_num(k) + beta(m,k)*sqrt(p_d*sum(eta_eq(m,setdiff(1:K,k),n).*beta_uc(m,setdiff(1:K,k))));
                        k_tilde_den(k) = k_tilde_den(k) + p_d*(beta(m,k)^2)*sum(eta_eq(m,setdiff(1:K,k),n).*(beta_uc(m,setdiff(1:K,k))).^2);
                        theta_tilde_num(k) = theta_tilde_num(k) + p_d*(beta(m,k)^2)*sum(eta_eq(m,setdiff(1:K,k),n).*(beta_uc(m,setdiff(1:K,k))).^2);
                        theta_tilde_den(k) = theta_tilde_den(k) + p_d*beta(m,k)*sum(sqrt(eta_eq(m,setdiff(1:K,k),n)).*(beta_uc(m,setdiff(1:K,k))));                      
                    end
                end
            end
        end
    else
        ues_not_rearranged = setdiff((1+K_mmW):K,ue_rearranged);
        for m = 1:M
            for k = 1:K
                if ismember(m,Serv{k})
                    if ((k<=K_mmW) && (sub6ConnectionState(k) == 1))
                        for n = 1:num_sc_sub6
    %                         term = (N_AP*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*beta_uc(m,(1+K_mmW):K)*user_sc_alloc((1+K_mmW):K,n)));
                            term = num_sc_sub6*N_AP*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                                              
                            if term > 0
                                eta_eq(m,k,n) = p_fac/term;
                            end
                            N_tilde_num(k) = N_tilde_num(k) + sqrt(p_d*eta_eq(m,k,n))*beta_uc(m,k);
                            N_tilde_den(k) = N_tilde_den(k) + p_d*eta_eq(m,k,n)*(beta_uc(m,k))^2;
                            beta_tilde_num(k) = beta_tilde_num(k) + p_d*eta_eq(m,k,n)*(beta_uc(m,k))^2;
                            beta_tilde_den(k) = beta_tilde_den(k) + sqrt(p_d*eta_eq(m,k,n))*beta_uc(m,k);
                            k_tilde_num(k) = k_tilde_num(k) + beta(m,k)*sqrt(p_d*sum(eta_eq(m,setdiff(1:K,k),n).*beta_uc(m,setdiff(1:K,k))));
                            k_tilde_den(k) = k_tilde_den(k) + p_d*(beta(m,k)^2)*sum(eta_eq(m,setdiff(1:K,k),n).*(beta_uc(m,setdiff(1:K,k))).^2);
                            theta_tilde_num(k) = theta_tilde_num(k) + p_d*(beta(m,k)^2)*sum(eta_eq(m,setdiff(1:K,k),n).*(beta_uc(m,setdiff(1:K,k))).^2);
                            theta_tilde_den(k) = theta_tilde_den(k) + p_d*beta(m,k)*sum(sqrt(eta_eq(m,setdiff(1:K,k),n)).*(beta_uc(m,setdiff(1:K,k))));                      
                        end
                    elseif (k>K_mmW)
    %                     eta_eq(m,k) = 1./(N_AP*(N_UE_mmW*p_fac*beta_uc(m,1:K_mmW)+N_UE_sub6*sum(beta_uc(m,2:K))));
                        for n = 1:num_sc_sub6
    %                         term = (N_AP*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*beta_uc(m,(1+K_mmW):K)*user_sc_alloc((1+K_mmW):K,n)));
                            term = num_sc_sub6*N_AP*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                        
                            if term > 0
                                if ismember(k,ue_rearranged) 
                                    eta_eq(m,k,n) = p_fac_rearrange/term;
                                else
                                    eta_eq(m,k,n) = 1/term;
                                end
                            end
                            N_tilde_num(k) = N_tilde_num(k) + sqrt(p_d*eta_eq(m,k,n))*beta_uc(m,k);
                            N_tilde_den(k) = N_tilde_den(k) + p_d*eta_eq(m,k,n)*(beta_uc(m,k))^2;
                            beta_tilde_num(k) = beta_tilde_num(k) + p_d*eta_eq(m,k,n)*(beta_uc(m,k))^2;
                            beta_tilde_den(k) = beta_tilde_den(k) + sqrt(p_d*eta_eq(m,k,n))*beta_uc(m,k);
                            k_tilde_num(k) = k_tilde_num(k) + beta(m,k)*sqrt(p_d*sum(eta_eq(m,setdiff(1:K,k),n).*beta_uc(m,setdiff(1:K,k))));
                            k_tilde_den(k) = k_tilde_den(k) + p_d*(beta(m,k)^2)*sum(eta_eq(m,setdiff(1:K,k),n).*(beta_uc(m,setdiff(1:K,k))).^2);
                            theta_tilde_num(k) = theta_tilde_num(k) + p_d*(beta(m,k)^2)*sum(eta_eq(m,setdiff(1:K,k),n).*(beta_uc(m,setdiff(1:K,k))).^2);
                            theta_tilde_den(k) = theta_tilde_den(k) + p_d*beta(m,k)*sum(sqrt(eta_eq(m,setdiff(1:K,k),n)).*(beta_uc(m,setdiff(1:K,k))));                      
                        end
                    end
                end
            end
        end
    end
    N_tilde = N.*((N_tilde_num.^2)./N_tilde_den);
    beta_tilde = beta_tilde_num./beta_tilde_den;
    k_tilde = zeros(K,1);
    theta_tilde = zeros(K,1);
    for k = 1:K
        if(k_tilde_den(k) > 0)
            k_tilde(k) = (k_tilde_num(k))^2/(k_tilde_den(k));
        end
        if(theta_tilde_den(k) > 0)
            theta_tilde(k) = (theta_tilde_num(k)^2)/(theta_tilde_den(k));
        end
    end  
    pdTrue = GeneralizedGamma(1.37, 0.98, 1.60);
end