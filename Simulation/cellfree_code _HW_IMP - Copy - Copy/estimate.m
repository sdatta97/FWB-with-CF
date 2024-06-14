function [R_G_HAT_concat,C_G_HAT_concat,R_G,Wd_0,Wd_1,WdC_0,WdC_1,R_G_HAT,C_G_HAT] = estimate(AP,UE,N,R,p_pilot,pol,alpha,kappa,alpha_ADC,DA)
tau_p=UE;
if pol==1
 C_YY = zeros(N,N,AP,UE);
 Wd_0=zeros(N,N,AP);
 Wd_1=zeros(N,N,AP);
%  if DA==1
  for ap=1:AP
        for ue=1:UE
            Wd_0(:,:,ap)= diag(diag(p_pilot(ap,ue)*R(:,:,ap,ue)))+ Wd_0(:,:,ap);            
        end
  end
%  end
 for ap=1:AP
  for ue=1:UE
   for ue1=1:UE
        if rem((ue-ue1),tau_p)==0          
          C_YY(:,:,ap,ue) = tau_p*p_pilot(ap,ue1)*R(:,:,ap,ue1) + C_YY(:,:,ap,ue);
        else         
          C_YY(:,:,ap,ue) =C_YY(:,:,ap,ue) ;
        end
   end
%          if DA==1
          C_YY(:,:,ap,ue) =eye(N) +C_YY(:,:,ap,ue)+ kappa^2*Wd_0(:,:,ap)+(1-alpha_ADC)/alpha_ADC*((1+kappa^2)*Wd_0(:,:,ap)+eye(N))  ;
%          else
%            C_YY(:,:,ap,ue) =eye(N) +C_YY(:,:,ap,ue);  
%          end
%          
         R_G_HAT(:,:,ap,ue) = tau_p*p_pilot(ap,ue)*R(:,:,ap,ue)*pinv(C_YY(:,:,ap,ue))*(R(:,:,ap,ue))';

         C_G_HAT(:,:,ap,ue) = R(:,:,ap,ue) -  R_G_HAT(:,:,ap,ue);
         R_G(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = R(:,:,ap,ue);
         R_G_HAT_concat(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = R_G_HAT(:,:,ap,ue);
         C_G_HAT_concat(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = C_G_HAT(:,:,ap,ue);
         
 end
 end
 WdC_0=zeros(N*AP,N*AP);
 WdC_1=zeros(N*AP,N*AP);
   if DA==1
        for ue=1:UE
            WdC_0= diag(diag(p_pilot(1,ue)*R_G(:,:,ue)))+ WdC_0;           
        end
   end

elseif pol==2
    C_YY_1 = zeros(2*N,2*N,AP,UE);
    C_YY_2= zeros(2*N,2*N,AP,UE);
    R_G_HAT_concat=zeros(2*N*AP,2*N*AP,UE,2);
    R_G_HAT_concat=zeros(2*N*AP,2*N*AP,UE,2);
    R_G=zeros(2*N*AP,2*N*AP,UE,2);
    Wd_0=zeros(2*N,2*N,AP);
    Wd_1=zeros(2*N,2*N,AP);
    
for ap=1:AP
 for ue=1:UE
    R_G_1(:,:,ap,ue) = [(1-alpha)*R(:,:,ap,ue) zeros(N,N) ;zeros(N,N) alpha*R(:,:,ap,ue)];
    R_G_2(:,:,ap,ue) = [alpha*R(:,:,ap,ue) zeros(N,N); zeros(N,N) (1-alpha)*R(:,:,ap,ue)];
 end
end
% if DA==1
 for ap=1:AP
        for ue=1:UE
            Wd_0(:,:,ap)= diag(diag(p_pilot(ap,ue)*R_G_1(:,:,ap,ue)))+ Wd_0(:,:,ap);
            Wd_1(:,:,ap)= diag(diag(p_pilot(ap,ue)*R_G_2(:,:,ap,ue)))+ Wd_1(:,:,ap);
        end
 end
% end


for ap=1:AP
 for ue=1:UE
    for ue1=1:UE
        if rem(abs(ue-ue1),tau_p)==0      
         C_YY_1(:,:,ap,ue) = tau_p*p_pilot(ap,ue)*( R_G_1(:,:,ap,ue1))  +C_YY_1(:,:,ap,ue) ;
         C_YY_2(:,:,ap,ue) = tau_p*p_pilot(ap,ue)*( R_G_2(:,:,ap,ue1))  + C_YY_2(:,:,ap,ue);
         
        else
         C_YY_1(:,:,ap,ue) = C_YY_1(:,:,ap,ue) ;
         C_YY_2(:,:,ap,ue) =C_YY_2(:,:,ap,ue);
        
        end
  
    end
%     if DA==1
    C_YY_1(:,:,ap,ue) =eye(2*N) +C_YY_1(:,:,ap,ue)+ kappa^2*Wd_0(:,:,ap)+(1-alpha_ADC)/alpha_ADC*((1+kappa^2)*Wd_0(:,:,ap)+eye(2*N)) ;
    C_YY_2(:,:,ap,ue) = eye(2*N) + C_YY_2(:,:,ap,ue) + kappa^2*Wd_1(:,:,ap)+(1-alpha_ADC)/alpha_ADC*((1+kappa^2)*Wd_1(:,:,ap)+eye(2*N)) ;
%     else
%     C_YY_1(:,:,ap,ue) =eye(2*N) +C_YY_1(:,:,ap,ue);
%     C_YY_2(:,:,ap,ue) = eye(2*N) + C_YY_2(:,:,ap,ue)  ;  
%     end
    R_G_HAT_1(:,:,ap,ue) = tau_p*p_pilot(ap,ue)*R_G_1(:,:,ap,ue)*pinv(C_YY_1(:,:,ap,ue))*(R_G_1(:,:,ap,ue))';
    R_G_HAT_2(:,:,ap,ue) = tau_p*p_pilot(ap,ue)*R_G_2(:,:,ap,ue)*pinv(C_YY_2(:,:,ap,ue))*(R_G_2(:,:,ap,ue))';
    C_G_HAT_1(:,:,ap,ue) = R_G_1(:,:,ap,ue) -  R_G_HAT_1(:,:,ap,ue);
    C_G_HAT_2(:,:,ap,ue) = R_G_2(:,:,ap,ue) -  R_G_HAT_2(:,:,ap,ue);
%     R_G_HAT(:,:,ap,ue) = [R_G_HAT_1(:,:,ap,ue) zeros(N,N);zeros(N,N) R_G_HAT_2(:,:,ap,ue)];
%     C_G_HAT(:,:,ap,ue) = [C_G_HAT_1(:,:,ap,ue) zeros(N,N);zeros(N,N) C_G_HAT_2(:,:,ap,ue)];
    R_G(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N,ue,1) = R_G_1(:,:,ap,ue);
    R_G(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N,ue,2) = R_G_2(:,:,ap,ue);
    R_G_HAT_concat(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N,ue,1) = R_G_HAT_1(:,:,ap,ue);
    C_G_HAT_concat(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N,ue,1) = C_G_HAT_1(:,:,ap,ue);
    R_G_HAT_concat(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N,ue,2) = R_G_HAT_2(:,:,ap,ue);
    C_G_HAT_concat(((ap-1)*2*N+1):ap*2*N,((ap-1)*2*N+1):ap*2*N,ue,2) = C_G_HAT_2(:,:,ap,ue);
    R_G_HAT(:,:,1,ap,ue)=R_G_HAT_1(:,:,ap,ue);
    R_G_HAT(:,:,2,ap,ue)=R_G_HAT_2(:,:,ap,ue);
    C_G_HAT(:,:,1,ap,ue)  =C_G_HAT_1(:,:,ap,ue);
    C_G_HAT(:,:,2,ap,ue)  =C_G_HAT_2(:,:,ap,ue);
 end
end


WdC_0=zeros(2*N*AP,2*N*AP);
WdC_1=zeros(2*N*AP,2*N*AP);

  if DA==1    
        for ue=1:UE
            WdC_0= diag(diag(p_pilot(1,ue)*R_G(:,:,ue,1)))+ WdC_0;
            WdC_1= diag(diag(p_pilot(1,ue)*R_G(:,:,ue,2)))+ WdC_1;
        end
  end
end
   