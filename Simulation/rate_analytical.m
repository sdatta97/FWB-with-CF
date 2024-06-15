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
      e_new(k) = (1/(M*N_AP))*trace(conj(Theta(:,:,k))*pinv(T_inv)); 
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
zeta = zeros(K,K);
for i=1:K
    for j=1:K
        zeta(i,j) = trace(conj(Theta(:,:,i))*T_prime(:,:,j))/((1+e_new(i))^2*f_prime(j));
    end
end
delta = (N_AP*M)*((e_new.^2)./f_prime);

rate_dl = zeros(K,1);
MUI = zeros(K,1);
for k = 1:K
    for q = 1:K
        if (q~=k)
            MUI(k) = MUI(k) + (1/(N_AP*M))*zeta(k,q);
        end
    end
    rate_dl(k) = TAU_FAC*log2(1+delta(k)/(MUI(k)+1));   
end
end