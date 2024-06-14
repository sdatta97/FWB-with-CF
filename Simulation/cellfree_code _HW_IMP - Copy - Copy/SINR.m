function [rate1,Num0_k,Num1_k,Den0_k,Den1_k]= SINR(UE,AP,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C)

for ue=1:UE
    Num0_k(ue)=0;
    Num1_k(ue)=0;
    Den0_k(ue)=0;
    Den1_k(ue)=0;
    for ap=1:AP
        Num0_k(ue)=Num0_k(ue) +b_mk(ap,ue)*c_mk_n0(ap,ue);
        Num1_k(ue)=Num1_k(ue) +b_mk(ap,ue)*c_mk_n1(ap,ue);
        for ue1=1:UE
          Den0_k(ue) =Den0_k(ue) + b_mk(ap,ue1)^2*c_mk_d0(ap,ue,ue1);
          Den1_k(ue) =Den1_k(ue) + b_mk(ap,ue1)^2*c_mk_d1(ap,ue,ue1);
        end
    end
    Num0_k(ue)=Num0_k(ue)^2;
    Num1_k(ue)=Num1_k(ue)^2;
    Den0_k(ue)=Den0_k(ue)+1;
    Den1_k(ue)=Den1_k(ue)+1;
    SINR0(ue)=Num0_k(ue)/(Den0_k(ue)); % =Num0_k(ue)^2;%
    SINR1(ue)=Num1_k(ue)/(Den1_k(ue)); % Num1_k(ue)^2;%
    rate(ue)= (1-2*UE/T_C)*(log2(1+SINR0(ue))+ log2(1+SINR1(ue)))- beta*(sqrt(2*SINR0(ue)./(1+SINR0(ue))) + sqrt(2*SINR1(ue)./(1+SINR1(ue))) );
end
rate1=sum(rate);