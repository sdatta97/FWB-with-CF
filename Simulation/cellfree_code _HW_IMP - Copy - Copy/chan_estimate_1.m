function [R_G_HAT,C_G_HAT,R_G_HAT_concat,C_G_HAT_concat] = chan_estimate_1 (AP,UE,N,R,p)
 g_hat = zeros(N,AP,UE);
   g_LOS = zeros(N,AP,UE);
   Q = zeros(N,N,AP,UE);
   Q1 = zeros(N,N,AP,UE);
   Q2 =zeros(N*AP,N*AP);
   Q22=zeros(N*AP,N*AP);
   R_G_HAT_concat=zeros(N*AP,N*AP,UE);
   C_G_HAT_concat=zeros(N*AP,N*AP,UE);
   y1=zeros(N*AP,UE);
   y=zeros(N,AP,UE);
   C_YY = zeros(N,N,AP,UE);
    C_YY1 = zeros(N*AP,N*AP,UE);

   tau_p=5;
for ap=1:AP
 for ue=1:UE
   for ue1=1:UE
        if rem(abs(ue-ue1),tau_p)==0
   
%           y(:,ap,ue)= tau_p*sqrt(p(ap,ue1))*( (G(:,ap,ue1))) +y(:,ap,ue);
    
          C_YY(:,:,ap,ue) = tau_p^2*p(ap,ue1)*R(:,:,ap,ue1) + C_YY(:,:,ap,ue);
        else
%            y(:,ap,ue)=y(:,ap,ue); 
           C_YY(:,:,ap,ue) =C_YY(:,:,ap,ue) ;
        end
   end
%           y(:,ap,ue)=sqrt(tau_p)*(randn(N,1)+1i*randn(N,1))*1/sqrt(2) +y(:,ap,ue);
          C_YY(:,:,ap,ue) =tau_p*eye(N) +C_YY(:,:,ap,ue) ;
    
%          g_hat(:,ap,ue) = tau_p*sqrt(p(ap,ue))*R(:,:,ap,ue)*pinv(C_YY(:,:,ap,ue))*y(:,ap,ue) ;
         R_G_HAT(:,:,ap,ue) = tau_p^2*p(ap,ue)*R(:,:,ap,ue)*pinv(C_YY(:,:,ap,ue))*(R(:,:,ap,ue))';
%          g_hat(:,ap,ue) = R_G_HAT(:,:,ap,ue)*(randn(N,1)+1i*randn(N,1))*1/sqrt(2);%
         C_G_HAT(:,:,ap,ue) = R(:,:,ap,ue) -  R_G_HAT(:,:,ap,ue);
%          mrc_combiner(:,ap,ue)=  g_hat(:,ap,ue)/norm( g_hat(:,ap,ue));%norm( g_hat(:,ap,ue));%sqrt(trace(R_G_HAT(:,:,ap,ue)));%
         R_G_concat(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = R(:,:,ap,ue);
         R_G_HAT_concat(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = R_G_HAT(:,:,ap,ue);
         C_G_HAT_concat(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = C_G_HAT(:,:,ap,ue);
         
 end
end