function rate_dl = rate_analytical(params,sub6ConnectionState)
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
        R_g((m-1)*N_AP+1:m*N_AP,(m-1)*N_AP+1:m*N_AP,k) = R_gNB(1:N_AP,1:N_AP,m,k);
    end
end
for k = 1:K
    if (k<=K_mmW)
        Theta(:,:,k) = p_d*p_fac*(R_g(:,:,k));
    else
        Theta(:,:,k) = p_d*(R_g(:,:,k));
    end
end
tol = 1e-1;
err = 1;
e_init = zeros(K,1);
e_new = zeros(K,1);

while err > tol 
   T_inv = (1/(M*N_AP))*eye(M*N_AP); 
   for k = 1:K  
       if ((k<=K_mmW) && (sub6ConnectionState(k)==0))
           continue;             
       else
           T_inv = T_inv + (1/(M*N_AP))*conj(Theta(:,:,k))/(1 + e_init(k)); 
       end
   end  
   for k = 1:K
       if ((k<=K_mmW) && (sub6ConnectionState(k)==0))
           continue;
       else
           e_new(k) = (1/(M*N_AP))*trace(conj(Theta(:,:,k))*pinv(T_inv)); 
       end
    end
    err  = norm(e_new - e_init)/norm(e_new);
    e_init = e_new;    
end

T = pinv(T_inv);
J = zeros(K,K);
V = zeros(K,K);
g = zeros(K,1);    
for k = 1:K
    if ((k<=K_mmW) && (sub6ConnectionState(k)==0))
        continue;
    else
        for q = 1:K
            if ((q<=K_mmW) && (sub6ConnectionState(q)==0))
                continue;
            else
                J(k,q)  = (1/(N_AP*M))*trace(conj(Theta(:,:,k))*T*conj(Theta(:,:,q))*T)/((N_AP*M)*(1+e_new(q))^2);
                V(k,q)  = (1/(N_AP*M))*trace(conj(Theta(:,:,k))*T*conj(Theta(:,:,q))*T);
            end
        end
        g(k)   = (1/(N_AP*M))*trace(conj(Theta(:,:,k))*T*T) ;
    end
end
f_prime  = pinv(eye(K)-J)*(g);
e_prime  = pinv(eye(K)-J)*(V);
for k = 1:K
    if ((k<=K_mmW) && (sub6ConnectionState(k)==0))
        continue;
    else
        T_prime(:,:,k) = T*conj(Theta(:,:,k))*T;
        for q = 1:K
            if ((q<=K_mmW) && (sub6ConnectionState(q)==0))
                continue;
            else
                T_prime(:,:,k) = T_prime(:,:,k) + T*((1/(N_AP*M))*conj(Theta(:,:,q))*e_prime(q,k)/(1+e_new(q))^2 + (1/(N_AP*M))*conj(Theta(:,:,q))*e_prime(q,k)/(1+e_new(q))^2)*T;
            end
        end
    end
end
zeta = zeros(K,K);
for k=1:K
    if ((k<=K_mmW) && (sub6ConnectionState(k)==0))
        continue;
    else
        for q=1:K
            if ((q<=K_mmW) && (sub6ConnectionState(q)==0))
                continue;
            else
                zeta(k,q) = trace(conj(Theta(:,:,k))*T_prime(:,:,q))/((1+e_new(k))^2*f_prime(q));
            end
        end
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