function [SIM_SE,V_sq,Delta_sq,V,Delta,add_INT] = SIM_SE_cal(G,p,combiner,K,M,N,pol,R,alpha,R_G,G_HAT,kappa,alpha_ADC,SIM_add_INT)
tau_c =200;
% tau_p = K;
if pol==2
  R_original =0;
else
    R_original = R;
end
if pol==2
V = zeros(2,2,K);
Delta = zeros(2,2,K);
SINR = zeros(K,1);
Delta_sq = zeros(2,2,K);
V_sq = zeros(2,2,K);
G1= G(:,1,:,:);
G2=G(:,2,:,:);
R_G1=R_G(:,:,1);
R_G2=R_G(:,:,2);
% G11=R_G1*(randn(2*N*M,K) +1i*randn(2*N*M,K))*(1/sqrt(2));
% G22=R_G2*(randn(2*N*M,K) +1i*randn(2*N*M,K))*(1/sqrt(2));
G11=reshape(G1,[2*N*M,K]);
G22=reshape(G2,[2*N*M,K]);% zeros(2*N*M,K);%
blk_1=zeros(2*M*N,2*M*N);
for ap=1:M
       blk_1(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N) = ones(2*N);
end
for ue=1:K
    G3(:,:,ue)=[G11(:,ue) G22(:,ue)];
    HW_INT(:,:,ue)= sqrtm(SIM_add_INT(:,:,ue))*(randn(2*M*N,2)+1i*randn(2*M*N,2))*1/sqrt(2);
end
% combiner =G_HAT;
for ue = 1:K
        V_sq(:,:,ue) =alpha_ADC*sqrt(p)*transpose(G3(:,:,ue))*conj(combiner(:,:,ue));
        V(:,:,ue) = alpha_ADC^2*(p*transpose(G3(:,:,ue))*conj(combiner(:,:,ue))*(transpose(G3(:,:,ue))*conj(combiner(:,:,ue)))') ;
        for ue1=1:K
            
                Delta(:,:,ue) =alpha_ADC^2*(p*transpose(G3(:,:,ue))*conj(combiner(:,:,ue1))*(transpose(G3(:,:,ue))*conj(combiner(:,:,ue1)))')+Delta(:,:,ue) ;%(p*combiner(:,:,AP,ue1)'*G(:,:,AP,ue)*G(:,:,AP,ue)'*combiner(:,:,AP,ue1)) +Delta(:,:,ue);
            
        end
    add_INT(:,:,ue) = (kappa^2*alpha_ADC + (1-alpha_ADC)*alpha_ADC)*(transpose(G3(:,:,ue))*SIM_add_INT(:,:,ue)*conj(G3(:,:,ue)));%HW_INT(:,:,ue))*(transpose(G3(:,:,ue))* HW_INT(:,:,ue))'
    Delta(:,:,ue) = Delta(:,:,ue)- V(:,:,ue) ;%*combiner(:,:,AP,ue)'*combiner(:,:,AP,ue);%
    SINR(ue) =(1-2*K/tau_c)*(log2(det(eye(2) + V(:,:,ue)*pinv(eye(2)+Delta(:,:,ue)+add_INT(:,:,ue)))));
end
elseif pol==1
  V = zeros(K,1);
   V_sq = zeros(K,1);
Delta = zeros(K,1);
Delta_sq = zeros(K,1);
SINR = zeros(K,1);
G_original=reshape(G,[N*M,K]);
G_error=G_original-G_HAT;
G2 = reshape(G,[N*M,K]);
% combiner=reshape(combiner,[N,M,K]);
% combiner =G_HAT;
for ue = 1:K

        V(ue) =alpha_ADC^2*p*abs(transpose(G2(:,ue))*conj(combiner(:,ue)))^2;%trace(conj(R_G_HAT_concat(:,:,ue))*conj(pinv(Q2)))/N;%p*abs(transpose(G2(:,ue))*conj(combiner(:,ue)))^2;%trace(conj(R_G_HAT_concat(:,:,ue))*conj(Q2));%
        V_sq(ue)=alpha_ADC*sqrt(p)*(transpose(G2(:,ue))*conj(combiner(:,ue)));

        for ue1=1:K
                                       
            Delta(ue) = alpha_ADC^2*p*abs(transpose(G2(:,ue))*conj (combiner(:,ue1)))^2 +Delta(ue);
        end

    Delta(ue) = Delta(ue) -V(ue);% 
    add_INT(ue) = (kappa^2*alpha_ADC + (1-alpha_ADC)*alpha_ADC)*(transpose(G2(:,ue))*SIM_add_INT(:,:,ue)*conj(G2(:,ue)));
    SINR(ue) =(1-K/tau_c)*log2(1 +(V(ue)/(1+Delta(ue))));
end
end
SIM_SE = sum(SINR);

