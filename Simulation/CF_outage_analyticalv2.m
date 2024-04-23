function p_out = CF_outage_analyticalv2(params,ue_idx)
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
    rmin = params.r_min(1);
    % %Communication bandwidth (Hz)
    % B = params.Band;
    B = scs_sub6;

    %Noise figure (in dB)
    noiseFigure = 7;

    %Compute noise power (in dBm)
    noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
    noisevar = db2pow(noiseVariancedBm(1));
    % noisevar = 1;
    % BETA = params.BETA;
    BETA = (params.BETA)*noisevar;
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
    N_tilde_num = zeros(K_mmW,1);
    N_tilde_den = zeros(K_mmW,1);
    beta_tilde_num = zeros(K_mmW,1);
    beta_tilde_den = zeros(K_mmW,1);
    k_tilde_num = zeros(K_mmW,1);
    k_tilde_den = zeros(K_mmW,1);
    theta_tilde_num = zeros(K_mmW,1);
    theta_tilde_den = zeros(K_mmW,1);    
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
                    end
                end
            end
        end
    end
    for k = 1:K_mmW
        if (sub6ConnectionState(k) == 1)
            for m = 1:M
                if ismember(m,Serv{k})
                    for n = 1:num_sc_sub6
                        N_tilde_num(k) = N_tilde_num(k) + sqrt(p_d*eta_eq(m,k,n))*beta_uc(m,k);
                        N_tilde_den(k) = N_tilde_den(k) + p_d*eta_eq(m,k,n)*(beta_uc(m,k))^2;
                        beta_tilde_num(k) = beta_tilde_num(k) + p_d*eta_eq(m,k,n)*(beta_uc(m,k))^2;
                        beta_tilde_den(k) = beta_tilde_den(k) + sqrt(p_d*eta_eq(m,k,n))*beta_uc(m,k);
                        k_tilde_num(k) = k_tilde_num(k) + BETA(m,k)*sqrt(p_d*sum(eta_eq(m,setdiff(1:K,k),n).*beta_uc(m,setdiff(1:K,k))));
                        k_tilde_den(k) = k_tilde_den(k) + p_d*(BETA(m,k)^2)*sum(eta_eq(m,setdiff(1:K,k),n).*(beta_uc(m,setdiff(1:K,k))).^2);
                        theta_tilde_num(k) = theta_tilde_num(k) + p_d*(BETA(m,k)^2)*sum(eta_eq(m,setdiff(1:K,k),n).*(beta_uc(m,setdiff(1:K,k))).^2);
                        theta_tilde_den(k) = theta_tilde_den(k) + p_d*BETA(m,k)*sum(sqrt(eta_eq(m,setdiff(1:K,k),n)).*(beta_uc(m,setdiff(1:K,k))));                      
                    end
                end
            end
        end
    end

    N_tilde = N.*((N_tilde_num.^2)./N_tilde_den);
    beta_tilde = beta_tilde_num./beta_tilde_den;
    k_tilde = zeros(K_mmW,1);
    theta_tilde = zeros(K_mmW,1);
    for k = 1:K_mmW
        if(k_tilde_den(k) > 0)
            k_tilde(k) = N*(k_tilde_num(k))^2/(k_tilde_den(k));
        end
        if(theta_tilde_den(k) > 0)
            theta_tilde(k) = theta_tilde_num(k)/theta_tilde_den(k);
        end
    end
    pdTrue = GeneralizedGamma(1/beta_tilde(ue_idx), 2, 0.5*N_tilde(ue_idx));
    if k_tilde(ue_idx)>0
        % y = (10.^(0:1:10)).*noisevar;
        mu = k_tilde(ue_idx)*theta_tilde(ue_idx);
        % sigma = sqrt(k_tilde(ue_idx))*theta_tilde(ue_idx);
        % y = (mu-sigma):noisevar:(mu+sigma);
        % y =  0:noisevar:(2*mu);
        y = 0:1e-10:0.001;
        stepsize = 1e-10;
        fy = gampdf(y,k_tilde(k),theta_tilde(k));
    else
        y = 0;
        stepsize = 1;
        fy = 1/noisevar;
    end
    x = rmin*(y+noisevar);
    % f = pdTrue.pdf(x);
    F = pdTrue.cdf(x);
    p_out = sum(fy.*F)*stepsize;
end