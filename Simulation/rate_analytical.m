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
R_g = zeros(M*N_AP,M*N_AP,K);
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
for k = 1:K
    for m = 1:M
        R_g((m-1)*N_AP+1:m*N_AP,(m-1)*N_AP+1:m*N_AP,k) = R_gNB(:,:,m,k);
    end
end
for k = 1:K
    Theta(:,:,k) = p_d*(R_g(:,:,k));
end
tol = 1e-1;
err = 1;
e_init = zeros(K,1);
e_new = zeros(K,1);

while err > tol 
   T_inv = (1/(M*N_AP))*eye(M*N_AP); 
   for k = 1:K         
       T_inv = T_inv + (1/(M*N_AP))*conj(Theta(:,:,k))/(1 + e_init(k));              
   end  
   for k = 1:K
      e_new(k) = (1/(M*N_AP))*trace(conj(Theta(:,:,k))*pinv(T_inv1)); 
    end
    err  = norm(e_new - e_init)/norm(e_new);
    e_init = e_new;    
end
T_inv_m = zeros(N_AP,N_AP,M);
T_inv_block=zeros(N_AP*M,N_AP*M);
for m=1:M
   T_inv_m(:,:,m) = T_inv(((m-1)*N_AP+1):m*N,((m-1)*N_AP+1):m*N_AP);
end
for m=1:M
   T_inv_block(((m-1)*N+1):m*N,((m-1)*N+1):m*N) = T_inv_m(:,:,m);
end

J = zeros(K,K);
V = zeros(K,K);
g = zeros(K,1);
for k = 1:K           
    Theta_k = p_d*R_g(:,:,k);
    Theta(:,:,k)=Theta_k;
end       
for m = 1:K
    for n = 1:K
      J(m,n)  = (1/(N_AP*M))*trace(conj(Theta(:,:,m))*pinv(T_inv)*conj(Theta(:,:,n))*pinv(T_inv))/((N_AP*M)*(1+e_new(n))^2);
      V(m,n)  = (1/(N_AP*M))*trace(conj(Theta(:,:,m))*pinv(T_inv)*conj(Theta(:,:,n))*pinv(T_inv));
    end
    g(m)   = (1/(N_AP*M))*trace(conj(Theta(:,:,m))*pinv(T_inv)*pinv(T_inv)) ;             
end
f_prime  = pinv(eye(K)-J)*(g);
e_prime  = pinv(eye(K)-J)*(V);
for k = 1:K
    T_prime(:,:,k) = pinv(T_inv)*conj(Theta(:,:,k))*pinv(T_inv);
    for q = 1:K
        T_prime(:,:,k) = T_prime(:,:,k) + pinv(T_inv)*((1/(N_AP*M))*conj(Theta(:,:,q))*e_prime(q,k)/(1+e(q))^2 + (1/(N_AP*M))*conj(Theta(:,:,q))*e_prime(q,k)/(1+e(q))^2)*pinv(T_inv);
    end
end
% T_prime_m = zeros(N_AP,N_AP,M);
% for k=1:K
%     for m=1:M
%        T_prime_m(:,:,m,k) = T_prime(((m-1)*N_AP+1):m*N_AP,((m-1)*N_AP+1):m*N_AP,k);
%     end
% end
% for k=1:K
%     M1(k) =(1/p_d)*((1/(N_AP*M))*f_prime(k)/((1+e_11(k))^2));
% end

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