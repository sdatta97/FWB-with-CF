
function [SINR_TH_UatF,M1 ,SIM_add_INT] = THH_UATF(p_data,p_opt,UE,AP,N ,R_G_HAT,pol,C_G_HAT,kappa,alpha_ADC,WDC_0,WDC_1,p_pilot,DA,beta,R_G_hat,C_G_hat)
T_C=200;
WDC_0=WDC_0*(p_data(1)/p_pilot(1,1));
WDC_1=WDC_1*(p_data(1)/p_pilot(1,1));
% beta=qfuncinv(decod_err_prob)*log2(exp(1))*sqrt(2)/sqrt(tau_b*B);
b_mk=ones(AP,UE);

if pol==2

SIGMA_N_BS=1;

R_G_HAT_1 =R_G_HAT(:,:,:,1);
R_G_HAT_2 =R_G_HAT(:,:,:,2);
C_G_HAT_1 =C_G_HAT(:,:,:,1);
C_G_HAT_2 =C_G_HAT(:,:,:,2);
% SINR_TH_UatF                  = zeros(1,K);
blk_1=zeros(2*AP*N,2*AP*N);
for ap=1:AP
       blk_1(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N) = ones(2*N);
end

tolerance = 10^(-7);
error2     = 1;
error1=1;


SIGMA_N_BS=1;

p_ue_sum_C_U1 = zeros(AP*2*N,AP*2*N);
p_ue_sum_C_U2 = zeros(AP*2*N,AP*2*N);
 for ue = 1:UE
   p_ue_sum_C_U1= p_ue_sum_C_U1 + p_data(ue)*C_G_HAT_1(:,:,ue);%
   p_ue_sum_C_U2= p_ue_sum_C_U2 + p_data(ue)*C_G_HAT_2(:,:,ue);
 end

tolerance = 10^(-9);

error1=1;
e_init2    = zeros(UE,1);
e_init1    = zeros(UE,1);
e_11       = zeros(UE,1);
e_22       = zeros(UE,1);

for ue = 1:UE
    Theta1(:,:,ue) = p_data(ue)*(R_G_HAT_1(:,:,ue));
    Theta2(:,:,ue) = p_data(ue)*(R_G_HAT_2(:,:,ue));
end

while max(error1,error2) > tolerance 
   if DA==1
    T_inv1 = (1/(2*AP*N))*conj( p_ue_sum_C_U1+p_ue_sum_C_U2+(SIGMA_N_BS^2)*eye(AP*2*N)*(1+(1-alpha_ADC)/alpha_ADC)+ (kappa^2+(kappa^2+1)*(1-alpha_ADC)/alpha_ADC)*(WDC_0+WDC_1)); 
   else
    T_inv1 = (1/(2*AP*N))*conj( (SIGMA_N_BS^2)*eye(AP*2*N));    
   end
    
     for ue = 1:UE
         T_inv1 = T_inv1 + (1/(2*AP*N))*conj(Theta1(:,:,ue))/(1 + e_init1(ue))+ (1/(2*AP*N))*conj(Theta2(:,:,ue))/(1 + e_init2(ue));
        
     end
%     
    for ue = 1:UE
      e_11(ue) = (1/(2*AP*N))*trace(conj(Theta1(:,:,ue))*pinv(T_inv1)); 
      e_22(ue) = (1/(2*AP*N))*trace(conj(Theta2(:,:,ue))*pinv(T_inv1)); 
      e(2*ue-1)=e_11(ue);
      e(2*ue)=e_22(ue);
    end
    error1  = norm(e_11 - e_init1)/norm(e_11);
    error2  = norm(e_22 - e_init2)/norm(e_22);
    e_init1 = e_11;
    e_init2 = e_22;
    
end
T_inv1_m = zeros(2*N,2*N,AP);
T_inv1_block=zeros(2*N*AP,2*N*AP);
for ap=1:AP
       T_inv1_m(:,:,ap) = T_inv1(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N);
end
for ap=1:AP
       T_inv1_block(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N) = T_inv1_m(:,:,ap);
end
for ap=1:AP
    for ue=1:UE
        e_111(ap,ue) = p_data(ue)*(1/(2*AP*N))*trace(conj(R_G_hat(:,:,1,ap,ue))*pinv(T_inv1_m(:,:,ap))); 
        e_222(ap,ue) = p_data(ue)*(1/(2*AP*N))*trace(conj(R_G_hat(:,:,2,ap,ue))*pinv(T_inv1_m(:,:,ap)));
    end
end

      J = zeros(UE,UE);
      V = zeros(UE,UE);
      v = zeros(UE,1);
      for m = 1:UE
%           
            Theta1_m = p_data(m)*(R_G_HAT_1(:,:,m));
            Theta2_m = p_data(m)*(R_G_HAT_2(:,:,m));
            Theta(:,:,2*m-1)=Theta1_m;
            Theta(:,:,2*m)=Theta2_m;
      end
%       
       for m = 1:2*UE
%           
%             Theta1_m = p*(R_G_HAT_1(:,:,m));
%             Theta2_m = p*(R_G_HAT_2(:,:,m));
          
            for n = 1:2*UE
%               
%               Theta1_n = p*(R_G_HAT_1(:,:,n));
%               Theta2_n = p*(R_G_HAT_2(:,:,n));
% %               
%                 J1(m,n)  = (1/(AP*N))*trace(conj(Theta1_m)*pinv(T_inv1)*conj(Theta1_n)*pinv(T_inv1))/((AP*N)*(1+e_11(n))^2);
%                 J2(m,n)  = (1/(AP*N))*trace(conj(Theta2_m)*pinv(T_inv2)*conj(Theta2_n)*pinv(T_inv2))/((AP*N)*(1+e_22(n))^2);
%                 V1(m,n)  = (1/(AP*N))*trace(conj(Theta1_m)*pinv(T_inv1)*conj(Theta1_n)*pinv(T_inv1));
%                 V2(m,n)  = (1/(AP*N))*trace(conj(Theta2_m)*pinv(T_inv2)*conj(Theta2_n)*pinv(T_inv2));
                  J(m,n)  = (1/(2*AP*N))*trace(conj(Theta(:,:,m))*pinv(T_inv1)*conj(Theta(:,:,n))*pinv(T_inv1))/((2*AP*N)*(1+e(n))^2);
                  V(m,n)  = (1/(2*AP*N))*trace(conj(Theta(:,:,m))*pinv(T_inv1)*conj(Theta(:,:,n))*pinv(T_inv1));
            end
%               
%               v1(m,1)  = (1/(AP*N))*trace(conj(Theta1_m)*inv(T_inv1)*inv(T_inv1))  ;
%               v2(m,1)  = (1/(AP*N))*trace(conj(Theta2_m)*inv(T_inv2)*inv(T_inv2)) ;
                v(m)   = (1/(2*AP*N))*trace(conj(Theta(:,:,m))*inv(T_inv1)*inv(T_inv1)) ;
             
       end
%      
%        e_prime1  = inv(eye(UE)-J1)*(v1);
%        e_prime2  = inv(eye(UE)-J2)*(v2);
         e_prime  = inv(eye(2*UE)-J)*(v);
%       
%        T_prime1     = inv(T_inv1)*(eye((AP*2*N)))*inv(T_inv1);
       
%        E_prime1           = inv(eye(UE)-J1)*(V1);
%        E_prime2           = inv(eye(UE)-J2)*(V2);
       E_prime           = inv(eye(2*UE)-J)*(V);
%      for ap=1:AP
         for r_1 = 1:UE
%               
%               T_prime3(:,:,r_1) = inv(T_inv1)*p*conj(R_G_HAT_1(:,:,r_1))*inv(T_inv1);
%               T_prime4(:,:,r_1) = inv(T_inv2)*p*conj(R_G_HAT_2(:,:,r_1))*inv(T_inv2);
                T_prime5(:,:,2*r_1-1) = inv(T_inv1)*conj(Theta(:,:,2*r_1-1))*inv(T_inv1);
                T_prime5(:,:,2*r_1) = inv(T_inv1)*conj(Theta(:,:,2*r_1))*inv(T_inv1);
%                 T_prime55(:,:,2*r_1-1,ap) = inv(T_inv1_m(:,:,ap))*conj(R_G_hat(:,:,1,ap,2*r_1-1))*inv(T_inv1_m(:,:,ap));
%                 T_prime55(:,:,2*r_1,ap) = inv(T_inv1_m(:,:,ap))*conj(R_G_hat(:,:,1,ap,2*r_1))*inv(T_inv1_m(:,:,ap));
            for s_1 = 1:UE
%             
%               T_prime3(:,:,r_1) = T_prime3(:,:,r_1) + (inv(T_inv1)*((1/(AP*N))*p*conj(R_G_HAT_1(:,:,s_1))*E_prime1(s_1,r_1)/(1+e_11(s_1))^2 + (1/(AP*N))*p*conj(R_G_HAT_2(:,:,s_1))*E_prime2(s_1,r_1)/(1+e_22(s_1))^2)*inv(T_inv1));
%               T_prime4(:,:,r_1) = T_prime4(:,:,r_1) + (inv(T_inv2)*((1/(AP*N))*p*conj(R_G_HAT_1(:,:,s_1))*E_prime1(s_1,r_1)/(1+e_11(s_1))^2 + (1/(AP*N))*p*conj(R_G_HAT_2(:,:,s_1))*E_prime2(s_1,r_1)/(1+e_22(s_1))^2)*inv(T_inv2));
              T_prime5(:,:,2*r_1-1) = T_prime5(:,:,2*r_1-1) + (inv(T_inv1)*((1/(2*AP*N))*conj(Theta(:,:,2*s_1-1))*E_prime(2*s_1-1,2*r_1-1)/(1+e(2*s_1-1))^2 + (1/(2*AP*N))*conj(Theta(:,:,2*s_1))*E_prime(2*s_1,2*r_1-1)/(1+e(2*s_1))^2)*inv(T_inv1));
              T_prime5(:,:,2*r_1) = T_prime5(:,:,2*r_1) + (inv(T_inv1)*((1/(2*AP*N))*conj(Theta(:,:,2*s_1-1))*E_prime(2*s_1-1,2*r_1)/(1+e(2*s_1-1))^2 + (1/(2*AP*N))*conj(Theta(:,:,2*s_1))*E_prime(2*s_1,2*r_1)/(1+e(2*s_1))^2)*inv(T_inv1));
%               T_prime55(:,:,2*r_1-1,ap) = T_prime55(:,:,2*r_1-1,ap) + (inv(T_inv1_m(:,:,ap))*((1/(2*AP*N))*conj(R_G_hat(:,:,1,ap,2*s_1-1))*E_prime(2*s_1-1,2*r_1-1)/(1+e(2*s_1-1))^2 + (1/(2*AP*N))*conj(R_G_hat(:,:,1,ap,2*s_1))*E_prime(2*s_1,2*r_1-1)/(1+e(2*s_1))^2)*inv(T_inv1_m(:,:,ap)));
%               T_prime55(:,:,2*r_1,ap) = T_prime55(:,:,2*r_1,ap) + (inv(T_inv1_m(:,:,ap))*((1/(2*AP*N))*conj(R_G_hat(:,:,1,ap,2*s_1-1))*E_prime(2*s_1-1,2*r_1)/(1+e(2*s_1-1))^2 + (1/(2*AP*N))*conj(R_G_hat(:,:,1,ap,2*s_1-1))*E_prime(2*s_1,2*r_1)/(1+e(2*s_1))^2)*inv(T_inv1_m(:,:,ap)));
            end
         end
%      end      
T_prime5_m = zeros(2*N,2*N,AP);
for r=1:UE
for ap=1:AP
       T_prime5_m(:,:,ap,2*r-1) = T_prime5(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N,2*r-1);
       T_prime5_m(:,:,ap,2*r) = T_prime5(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N,2*r);
end
end
for ue=1:UE
    M1(ue,1) =(1/p_data(ue))*((1/(2*AP*N))*e_prime(2*ue-1)/((1+e_11(ue))^2));
    M1(ue,2) =(1/p_data(ue))*((1/(2*AP*N))*e_prime(2*ue)/((1+e_22(ue))^2));
end
add_INT1=zeros(UE,1);
add_INT2=zeros(UE,1);
add_INT=zeros(2,2,UE);
for ue=1:UE
    for ue1=1:UE
    add_INT1(ue)=(1/(2*AP*N))*trace(conj(R_G_HAT_1(:,:,ue)+C_G_HAT_1(:,:,ue))*diag(diag(T_prime5(:,:,2*ue1-1))))/e_prime(2*ue1-1)+(1/(2*AP*N))*trace(conj(R_G_HAT_1(:,:,ue)+C_G_HAT_1(:,:,ue))*diag(diag(T_prime5(:,:,2*ue1))))/e_prime(2*ue1)+add_INT1(ue);
    add_INT2(ue)=(1/(2*AP*N))*trace(conj(R_G_HAT_2(:,:,ue)+C_G_HAT_2(:,:,ue))*diag(diag(T_prime5(:,:,2*ue1))))/e_prime(2*ue1)+(1/(2*AP*N))*trace(conj(R_G_HAT_2(:,:,ue)+C_G_HAT_2(:,:,ue))*diag(diag(T_prime5(:,:,2*ue1-1))))/e_prime(2*ue1-1)+add_INT2(ue);
    
%     add_INT1(ue)= trace(conj(R_G_HAT_1(:,:,ue)+C_G_HAT_1(:,:,ue))*diag(diag(conj(R_G_HAT_1(:,:,ue1)))))/trace(R_G_HAT_1(:,:,ue1))+add_INT1(ue);
%     add_INT2(ue)= trace(conj(R_G_HAT_2(:,:,ue)+C_G_HAT_2(:,:,ue))*diag(diag(conj(R_G_HAT_2(:,:,ue1)))))/trace(R_G_HAT_2(:,:,ue1))+add_INT2(ue);
    end
    add_INT(:,:,ue)=(kappa^2*alpha_ADC+alpha_ADC*(1-alpha_ADC))*p_data(ue)*[add_INT1(ue) 0; 0 add_INT2(ue)];
end
SIM_add_INT = zeros(AP*2*N,AP*2*N,UE);
for ue=1:UE
    for ue1=1:UE
    SIM_add_INT(:,:,ue)=(1/(2*AP*N))*p_data(ue1)*diag(diag(T_prime5(:,:,2*ue1-1)))/e_prime(2*ue1-1)+(1/(2*AP*N))*p_data(ue1)*diag(diag(T_prime5(:,:,2*ue1)))/e_prime(2*ue1)+SIM_add_INT(:,:,ue);
%     SIM_add_INT(:,:,ue)= diag(diag(conj(R_G_HAT_1(:,:,ue1))))/trace(R_G_HAT_1(:,:,ue1))+ diag(diag(conj(R_G_HAT_2(:,:,ue1))))/trace(R_G_HAT_2(:,:,ue1))+SIM_add_INT(:,:,ue);
    end
end
des_sig1=zeros(UE,1);
des_sig2=zeros(UE,1);
mu_user_intf1=zeros(UE,UE);
mu_user_intf2=zeros(UE,UE);
BU_1= zeros(UE,1);
BU_2= zeros(UE,1);
for ue=1:UE
   
   
 for ap=1:AP
     des_sig1(ue) =b_mk(ap,ue)*alpha_ADC*(1/p_data(ue))*(e_111(ap,ue)/(1+e_11(ue)))/sqrt()+des_sig1(ue);%alpha_ADC^2*(1/p_data(ue))^2*(e_11(ue)/(1+e_11(ue)))^2/M1(ue,1);
     des_sig2(ue) =b_mk(ap,ue)*alpha_ADC*(1/p_data(ue))*(e_222(ap,ue)/(1+e_22(ue)))+des_sig2(ue);%alpha_ADC^2*(1/p_data(ue))^2*(e_22(ue)/(1+e_22(ue)))^2/M1(ue,2);
 end
   des_sig1(ue)=des_sig1(ue)^2/M1(ue,1);
   des_sig2(ue)=des_sig2(ue)^2/M1(ue,2);
%      des_sig(:,:,ue) = alpha_ADC^2*[des_sig1(ue) 0; 0 des_sig2(ue)];
 for ap=1:AP
     for ue1=1:UE 
         if ue1~=ue
%            mu_user_intf1(ue,ue1)= (1/p_data(ue)^2)*(alpha_ADC^2*((1/(2*AP*N))^2*trace(conj(Theta1(:,:,ue))*T_prime5(:,:,2*ue1-1))/((1+e_11(ue))^2*(1+e_11(ue1))^2)/M1(ue1,1) + (1/(2*AP*N))^2*trace(p_data(ue)*conj(C_G_HAT_1(:,:,ue))*T_prime5(:,:,2*ue1-1))/((1+e_11(ue1))^2)/M1(ue1,1) + (1/(2*AP*N))^2*trace(conj(Theta1(:,:,ue))*T_prime5(:,:,2*ue1))/((1+e_11(ue))^2*(1+e_22(ue1))^2)/M1(ue1,2) + (1/(2*AP*N))^2*trace(p_data(ue)*conj(C_G_HAT_1(:,:,ue))*T_prime5(:,:,2*ue1))/((1+e_22(ue1))^2)/M1(ue1,2)));
%            mu_user_intf2(ue,ue1)=(1/p_data(ue)^2)* (alpha_ADC^2*((1/(2*AP*N))^2*trace(conj(Theta2(:,:,ue))*T_prime5(:,:,2*ue1))/((1+e_22(ue))^2*(1+e_22(ue1))^2)/M1(ue1,2) + (1/(2*AP*N))^2*trace(p_data(ue)*conj(C_G_HAT_2(:,:,ue))*T_prime5(:,:,2*ue1))/((1+e_22(ue1))^2)/M1(ue1,2)+ (1/(2*AP*N))^2*trace(conj(Theta2(:,:,ue))*T_prime5(:,:,2*ue1-1))/((1+e_22(ue))^2*(1+e_11(ue1))^2)/M1(ue1,1) + (1/(2*AP*N))^2*trace(p_data(ue)*conj(C_G_HAT_2(:,:,ue))*T_prime5(:,:,2*ue1-1))/((1+e_11(ue1))^2)/M1(ue1,1)));
           mu_user_intf1(ue,ue1)=mu_user_intf1(ue,ue1)+ (1/p_data(ue))*b_mk(ap,ue1)*(alpha_ADC^2*((1/(2*AP*N))^2*trace(conj(R_G_hat(:,:,1,ap,ue))*T_prime5_m(:,:,ap,2*ue1-1))/((1+e_11(ue))^2*(1+e_11(ue1))^2)/M1(ue1,1) + (1/(2*AP*N))^2*trace(conj(C_G_hat(:,:,1,ap,ue))*T_prime5_m(:,:,ap,2*ue1-1))/((1+e_11(ue1))^2)/M1(ue1,1) + (1/(2*AP*N))^2*trace(conj(R_G_hat(:,:,1,ap,ue))*T_prime5_m(:,:,ap,2*ue1))/((1+e_11(ue))^2*(1+e_22(ue1))^2)/M1(ue1,2) + (1/(2*AP*N))^2*trace(conj(C_G_hat(:,:,1,ap,ue))*T_prime5_m(:,:,ap,2*ue1))/((1+e_22(ue1))^2)/M1(ue1,2)));
           mu_user_intf2(ue,ue1)=mu_user_intf2(ue,ue1)+(1/p_data(ue))* b_mk(ap,ue1)*(alpha_ADC^2*((1/(2*AP*N))^2*trace(conj(R_G_hat(:,:,2,ap,ue))*T_prime5_m(:,:,ap,2*ue1))/((1+e_22(ue))^2*(1+e_22(ue1))^2)/M1(ue1,2) + (1/(2*AP*N))^2*trace(conj(C_G_hat(:,:,2,ap,ue))*T_prime5_m(:,:,ap,2*ue1))/((1+e_22(ue1))^2)/M1(ue1,2)+ (1/(2*AP*N))^2*trace(conj(R_G_hat(:,:,2,ap,ue))*T_prime5_m(:,:,ap,2*ue1-1))/((1+e_22(ue))^2*(1+e_11(ue1))^2)/M1(ue1,1) + (1/(2*AP*N))^2*trace(conj(C_G_hat(:,:,2,ap,ue))*T_prime5_m(:,:,ap,2*ue1-1))/((1+e_11(ue1))^2)/M1(ue1,1)));
%          else
%            mu_user_intf1(ue)=mu_user_intf1(ue);
%            mu_user_intf2(ue)=mu_user_intf2(ue);
%          
         end
     end
     BU_1(ue) =(1/p_data(ue))*b_mk(ap,ue)*(alpha_ADC^2*(1/(2*AP*N))^2*trace(conj(C_G_hat(:,:,1,ap,ue))*T_prime5_m(:,:,ap,2*ue-1))/((1+e_11(ue))^2)/M1(ue,1)+  (1/(2*AP*N))^2*trace(conj(C_G_hat(:,:,1,ap,ue))*T_prime5_m(:,:,ap,2*ue))/((1+e_22(ue))^2)/M1(ue,2) + (1/(2*AP*N))^2*trace(conj(R_G_hat(:,:,1,ap,ue))*T_prime5_m(:,:,ap,2*ue))/((1+e_11(ue))^2*(1+e_22(ue))^2)/M1(ue,2));% +(1/(2*AP*N))^2*trace(conj(Theta1(:,:,ue))*T_prime5(:,:,2*ue-1))/((1+e_11(ue))^2)/M1(ue,1)-p*des_sig1(ue);
     BU_2(ue) =(1/p_data(ue))*b_mk(ap,ue)*(alpha_ADC^2*(1/(2*AP*N))^2*trace(conj(C_G_hat(:,:,2,ap,ue))*T_prime5_m(:,:,ap,2*ue))/((1+e_22(ue))^2)/M1(ue,2)+  (1/(2*AP*N))^2*trace(p_data(ue)*conj(C_G_hat(:,:,2,ap,ue))*T_prime5_m(:,:,ap,2*ue-1))/((1+e_11(ue))^2)/M1(ue,1) + (1/(2*AP*N))^2*trace(conj(R_G_hat(:,:,2,ap,ue))*T_prime5_m(:,:,ap,2*ue-1))/((1+e_11(ue))^2*(1+e_22(ue))^2)/M1(ue,1));% +(1/(2*AP*N))^2*trace(conj(Theta2(:,:,ue))*T_prime5(:,:,2*ue))/((1+e_22(ue))^2)/M1(ue,2)-p*des_sig2(ue);
%      BU(:,:,ue) = [ (1/p_data(ue)^2)*BU_1(ue) 0 ; 0 (1/p_data(ue)^2)*BU_2(ue)];
%      mu_user_intf(:,:,ue) = alpha_ADC^2*[ (1/p_data^2)*mu_user_intf1(ue) 0 ; 0 (1/p_data^2)*mu_user_intf2(ue)] ;
%      SINR_TH_UatF1(:,:,ue) = p_opt(ue)*des_sig(:,:,ue)*pinv( mu_user_intf(:,:,ue)+BU(:,:,ue)+eye(2)+add_INT(:,:,ue));
%      Rate(ue)=(1-2*UE/T_C)*log2(det(eye(2)+SINR_TH_UatF1(:,:,ue))) - qfuncinv(decod_err_prob)*log2(exp(1))*sqrt(2*SINR_TH_UatF1(1,1,ue)./(1+SINR_TH_UatF1(1,1,ue))./(tau_b*B))- qfuncinv(decod_err_prob)*log2(exp(1))*sqrt(2*SINR_TH_UatF1(2,2,ue)./(1+SINR_TH_UatF1(2,2,ue))./(tau_b*B));
 end
 end
Den1= zeros(UE,1);
Den2=zeros(UE,1);
for ue=1:UE
    Num1(ue)=p_opt(ue)*des_sig1(ue);
    Num2(ue)=p_opt(ue)*des_sig2(ue);
    for ue1=1:UE
        if ue1~=ue
            Den1(ue)= p_opt(ue1)*mu_user_intf1(ue,ue1)+Den1(ue);
            Den2(ue)= p_opt(ue1)*mu_user_intf2(ue,ue1)+Den2(ue);
        end
    end
    Den1(ue) = 1+ Den1(ue)+ p_opt(ue)*BU_1(ue)+add_INT(1,1,ue);
    Den2(ue) = 1+ Den2(ue)+ p_opt(ue)*BU_2(ue)+add_INT(2,2,ue);
    Num(:,:,ue)= [Num1(ue) 0; 0 Num2(ue)];
    Den(:,:,ue)= [Den1(ue) 0; 0 Den2(ue)];
    SINR_TH_UatF1(:,:,ue) = Num(:,:,ue)*pinv( Den(:,:,ue));
    Rate(ue)=(1-2*UE/T_C)*log2(det(eye(2)+SINR_TH_UatF1(:,:,ue))) - beta*(sqrt(2*SINR_TH_UatF1(1,1,ue)./(1+SINR_TH_UatF1(1,1,ue))) + sqrt(2*SINR_TH_UatF1(2,2,ue)./(1+SINR_TH_UatF1(2,2,ue))) );
end
SINR_TH_UatF = sum(Rate); 
else

 SIGMA_N_BS=1;
p_ue_sum_C_U = zeros(N,N,AP);
p_ue_sum_C_U4 = zeros(AP*N,AP*N);

 for ue = 1:UE
   p_ue_sum_C_U4= p_ue_sum_C_U4 + p_data*C_G_HAT(:,:,ue);%
 end
tolerance = 10^(-7);

error=1;
e_init    = zeros(UE,AP);
e_init1    = zeros(UE,1);
e_1       = zeros(UE,AP);
e_11       = zeros(UE,1);


for ue = 1:UE
    Theta(:,:,ue) = p_data*(R_G_HAT(:,:,ue));
   
end

while max((error)) > tolerance 
   if DA==1
    T_inv4 = (1/(AP*N))*conj( p_ue_sum_C_U4+(SIGMA_N_BS^2)*eye(AP*N)*(1+(1-alpha_ADC)/alpha_ADC)+ (kappa^2+(kappa^2+1)*(1-alpha_ADC)/alpha_ADC)*(WDC_0) );
   else
    T_inv4 = (1/(AP*N))*conj((SIGMA_N_BS^2)*eye(AP*N));
   end
    
     for ue = 1:UE
         T_inv4 = T_inv4 + (1/(AP*N))*conj(Theta(:,:,ue))/(1 + e_init1(ue));
     end
%     
    for ue = 1:UE
      e_11(ue) = (1/(AP*N))*trace(conj(Theta(:,:,ue))*pinv(T_inv4)); 
      
    end
    error  = norm(e_11 - e_init1)/norm(e_11);
    e_init1 = e_11;
end

      J = zeros(UE,UE);
      V = zeros(UE,UE);
      v = zeros(UE,AP);
%       
       for m = 1:UE
%           
            Theta_m = p_data*(R_G_HAT(:,:,m));
          
            for n = 1:UE
%               
              Theta_n = p_data*(R_G_HAT(:,:,n));
              
%               
                J1(m,n)  = (1/(AP*N))*trace(conj(Theta_m)*pinv(T_inv4)*conj(Theta_n)*pinv(T_inv4))/((AP*N)*(1+e_11(n))^2);
                V1(m,n)  = (1/(AP*N))*trace(conj(Theta_m)*pinv(T_inv4)*conj(Theta_n)*pinv(T_inv4));
            end
%               
              v1(m,1)  = (1/(AP*N))*trace(conj(Theta_m)*inv(T_inv4)*inv(T_inv4));
             
       end
%      
       e_prime1  = inv(eye(UE)-J1)*v1;
%        e_prime1(m) = (1/N)*trace(Theta_m*T_prime1);
%       
       T_prime1     = inv(T_inv4)*(eye((AP*N)))*inv(T_inv4);
%        
       E_prime1           = inv(eye(UE)-J1)*V1;
      
%       
       T_prime3           = zeros(N*AP,N*AP,UE);
%     
         for r_1 = 1:UE
%               
              T_prime3(:,:,r_1) = inv(T_inv4)*p_data*conj(R_G_HAT(:,:,r_1))*inv(T_inv4);
            for s_1 = 1:UE
%             
              T_prime3(:,:,r_1) = T_prime3(:,:,r_1) + (inv(T_inv4)*((1/(AP*N))*p_data*conj(R_G_HAT(:,:,s_1))*E_prime1(s_1,r_1)/(1+e_11(s_1))^2)*inv(T_inv4));
            end
         end
%        

for ue=1:UE
    M1(ue) =(1/p_data)*((1/(AP*N))*e_prime1(ue)/((1+e_11(ue))^2));
end
add_INT=zeros(UE,1);
for ue=1:UE
    for ue1=1:UE
    add_INT(ue)=(kappa^2*alpha_ADC+alpha_ADC*(1-alpha_ADC))*p_data*(1/(AP*N))*trace(conj(R_G_HAT(:,:,ue)+C_G_HAT(:,:,ue))*diag(diag(T_prime3(:,:,ue1))))/e_prime1(ue1)+add_INT(ue);    
% %     add_INT(ue)= trace(conj(R_G_HAT(:,:,ue)+C_G_HAT(:,:,ue))*diag(diag(conj(R_G_HAT(:,:,ue1)))))/trace(R_G_HAT(:,:,ue1))+add_INT(ue);
    end
%     add_INT(ue)=(kappa^2*alpha_ADC+alpha_ADC*(1-alpha_ADC))*p_data*(1/(AP*N))*trace(conj(R_G_HAT(:,:,ue)+C_G_HAT(:,:,ue))*diag(diag(T_prime3(:,:,ue))))/e_prime1(ue);
end
SIM_add_INT = zeros(AP*N,AP*N,UE);
for ue=1:UE
    for ue1=1:UE
    SIM_add_INT(:,:,ue)=(1/(AP*N))*p_data*diag(diag(T_prime3(:,:,ue1)))/e_prime1(ue1)+SIM_add_INT(:,:,ue);
%     SIM_add_INT(:,:,ue)= diag(diag(conj(R_G_HAT(:,:,ue1))))/trace(R_G_HAT(:,:,ue1))+SIM_add_INT(:,:,ue);
    end
end
for ue=1:UE
%     
    mu_user_intf(ue) = 0;

     des_sig(ue) =alpha_ADC^2*(1/p_data)*(e_11(ue)/(1+e_11(ue)))^2/M1(ue);
     for ue1=1:UE 
         if ue1~=ue
%         
           mu_user_intf(ue)= mu_user_intf(ue)+(1/(AP*N))^2*trace(conj(Theta(:,:,ue))*T_prime3(:,:,ue1))/((1+e_11(ue))^2*(1+e_11(ue1))^2)/M1(ue1)+(1/(AP*N))^2*trace(p_data*conj(C_G_HAT(:,:,ue))*T_prime3(:,:,ue1))/((1+e_11(ue1))^2)/M1(ue1);
         else
        
           mu_user_intf(ue)=mu_user_intf(ue);
%          
         end
     end
     
    BU(ue) =alpha_ADC^2*(1/p_data)*(1/(AP*N))^2*trace(p_data*conj(C_G_HAT(:,:,ue))*T_prime3(:,:,ue))/((1+e_11(ue))^2)/M1(ue);

 
 mu_user_intf(ue) = alpha_ADC^2*(1/p_data)*mu_user_intf(ue);
 SINR_TH_UatF1(ue) =(1-UE/T_C)*log2(1+add_INT(ue));%(1-UE/T_C)*log2(1+des_sig(ue)/(1+mu_user_intf(ue)+BU(ue)+add_INT(ue)));
end
SINR_TH_UatF = sum(SINR_TH_UatF1);   
end 
end























