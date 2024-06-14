function [g_hat,mmse_combiner,mrc_combiner,G_HAT, H] =chan_estimate (AP,UE,N,G,R,alpha,p_pilot,pol,p_transmit,M,R_G_HAT_concat,C_G_HAT_concat,a1,Norm_factor_true,kappa,alpha_ADC,Wd_0,Wd_1,WdC_0,WdC_1,DA)
tau_p=UE/1;
sigma=1;
if pol==2
g_hat = zeros(2*N,2,AP,UE);
G1=G(:,1,:,:);
G2=G(:,2,:,:);
H1=reshape(G1,[2*N*AP,UE]);
H2=reshape(G2,[2*N*AP,UE]);
for ue=1:UE
    H(:,:,ue)=[H1(:,ue) H2(:,ue)];
end
y1=zeros(2*N,AP,UE);
y2=zeros(2*N,AP,UE);
C_YY_1 = zeros(2*N,2*N,AP,UE);
C_YY_2 =zeros(2*N,2*N,AP,UE);
Q_1 = zeros(2*N*AP,2*N*AP);
Q_2 =zeros(2*N*AP,2*N*AP);


for ap=1:AP
 for ue=1:UE
    R_G_1(:,:,ap,ue) = [(1-alpha)*R(:,:,ap,ue) zeros(N,N) ;zeros(N,N) alpha*R(:,:,ap,ue)];
    R_G_2(:,:,ap,ue) = [alpha*R(:,:,ap,ue) zeros(N,N); zeros(N,N) (1-alpha)*R(:,:,ap,ue)];
 end
end
for ap=1:AP
 for ue=1:UE
    for ue1=1:UE
        if rem(abs(ue-ue1),tau_p)==0
         y1(:,ap,ue)= tau_p*sqrt(p_pilot(ap,ue))*( (G(:,1,ap,ue1))) +y1(:,ap,ue);
         y2(:,ap,ue) = tau_p*sqrt(p_pilot(ap,ue))*( (G(:,2,ap,ue1))) + y2(:,ap,ue);
         C_YY_1(:,:,ap,ue) = tau_p^1*p_pilot(ap,ue)*( R_G_1(:,:,ap,ue1))  +C_YY_1(:,:,ap,ue) ;
         C_YY_2(:,:,ap,ue) = tau_p^1*p_pilot(ap,ue)*( R_G_2(:,:,ap,ue1))  + C_YY_2(:,:,ap,ue);
         
        else
          y1(:,ap,ue)= y1(:,ap,ue);
         y2(:,ap,ue) =  y2(:,ap,ue);
         C_YY_1(:,:,ap,ue) = C_YY_1(:,:,ap,ue) ;
         C_YY_2(:,:,ap,ue) =C_YY_2(:,:,ap,ue);
        
        end
  
    end
    if DA==1
    y1_RF(:,ap,ue) = sqrtm(alpha_ADC^2*kappa^2*tau_p*Wd_0(:,:,ap))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2);
    y1_ADC(:,ap,ue) = sqrtm(tau_p*alpha_ADC*(1-alpha_ADC)*((1+kappa^2)*Wd_0(:,:,ap)+sigma^2*eye(2*N)))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2);
    y2_RF(:,ap,ue) = sqrtm(alpha_ADC^2*kappa^2*tau_p*Wd_1(:,:,ap))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2);
    y2_ADC(:,ap,ue) = sqrtm(tau_p*alpha_ADC*(1-alpha_ADC)*((1+kappa^2)*Wd_1(:,:,ap)+sigma^2*eye(2*N)))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2);
    y1(:,ap,ue)= sqrtm(alpha_ADC^2*tau_p*eye(2*N))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2) +alpha_ADC*y1(:,ap,ue) +y1_RF(:,ap,ue)+y1_ADC(:,ap,ue);
    y2(:,ap,ue) = sqrtm(alpha_ADC^2*tau_p*eye(2*N))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2) + alpha_ADC*y2(:,ap,ue)+y2_RF(:,ap,ue)+y2_ADC(:,ap,ue);
    else
    y1(:,ap,ue)= sqrtm(alpha_ADC^2*tau_p*eye(2*N))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2) +alpha_ADC*y1(:,ap,ue) ;
    y2(:,ap,ue) = sqrtm(alpha_ADC^2*tau_p*eye(2*N))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2) + alpha_ADC*y2(:,ap,ue);
    end
    C_YY_1(:,:,ap,ue) =eye(2*N) +C_YY_1(:,:,ap,ue)+kappa^2*Wd_0(:,:,ap)+(1-alpha_ADC)/alpha_ADC*((1+kappa^2)*Wd_0(:,:,ap)+eye(2*N)) ;
    C_YY_2(:,:,ap,ue) = eye(2*N) + C_YY_2(:,:,ap,ue) + kappa^2*Wd_1(:,:,ap)+(1-alpha_ADC)/alpha_ADC*((1+kappa^2)*Wd_1(:,:,ap)+eye(2*N));
%     else
%       y1(:,ap,ue)= sqrtm(alpha_ADC^2*tau_p*eye(2*N))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2) +alpha_ADC*y1(:,ap,ue);
%       y2(:,ap,ue) = sqrtm(alpha_ADC^2*tau_p*eye(2*N))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2) + alpha_ADC*y2(:,ap,ue);
%       C_YY_1(:,:,ap,ue) =eye(2*N) +C_YY_1(:,:,ap,ue);
%       C_YY_2(:,:,ap,ue) = eye(2*N) + C_YY_2(:,:,ap,ue);
%     end
    R_G_HAT_1(:,:,ap,ue) = tau_p^1*p_pilot(ap,ue)*R_G_1(:,:,ap,ue)*pinv(C_YY_1(:,:,ap,ue))*(R_G_1(:,:,ap,ue))';
    R_G_HAT_2(:,:,ap,ue) = tau_p^1*p_pilot(ap,ue)*R_G_2(:,:,ap,ue)*pinv(C_YY_2(:,:,ap,ue))*(R_G_2(:,:,ap,ue))';
    C_G_HAT_1(:,:,ap,ue) = R_G_1(:,:,ap,ue) -  R_G_HAT_1(:,:,ap,ue);
    C_G_HAT_2(:,:,ap,ue) = R_G_2(:,:,ap,ue) -  R_G_HAT_2(:,:,ap,ue);

    g_hat_1(:,ap,ue) = (1/alpha_ADC)*sqrt(p_pilot(ap,ue))*R_G_1(:,:,ap,ue)*pinv(C_YY_1(:,:,ap,ue))*y1(:,ap,ue);
    g_hat_2(:,ap,ue) = (1/alpha_ADC)*sqrt(p_pilot(ap,ue))*R_G_2(:,:,ap,ue)*pinv(C_YY_2(:,:,ap,ue))*y2(:,ap,ue); 
    mrc_combiner1(:,ap,ue)= g_hat_1(:,ap,ue)/sqrt(trace(R_G_1(:,:,ap,ue)));%norm(g_hat_1(:,ap,ue));%sqrt(trace(R_G_1(:,:,ap,ue)));
    mrc_combiner2(:,ap,ue)= g_hat_2(:,ap,ue)/sqrt(trace(R_G_2(:,:,ap,ue)));%norm(g_hat_2(:,ap,ue));%sqrt(trace(R_G_2(:,:,ap,ue)));
    g_hat(:,:,ap,ue)=[(g_hat_1(:,ap,ue) ) (g_hat_2(:,ap,ue) )];
%     mrc_combiner(:,:,ap,ue)=[mrc_combiner1(:,ap,ue)  mrc_combiner2(:,ap,ue) ];

  
 end
end
R_G_HAT(:,:,:,:,1)= R_G_HAT_1;
R_G_HAT(:,:,:,:,2)= R_G_HAT_2;
G_HAT1=reshape(g_hat_1,[2*N*AP,UE]);
G_HAT2=reshape(g_hat_2,[2*N*AP,UE]);

 for ue=1:UE
      
         Q_1 = p_transmit(1,ue)*( G_HAT1(:,ue)*G_HAT1(:,ue)' +C_G_HAT_concat(:,:,ue,1)+G_HAT2(:,ue)*G_HAT2(:,ue)'+C_G_HAT_concat(:,:,ue,2))  + Q_1;
         Q_2 = p_transmit(1,ue)*( G_HAT2(:,ue)*G_HAT2(:,ue)'+C_G_HAT_concat(:,:,ue,2)+G_HAT1(:,ue)*G_HAT1(:,ue)' +C_G_HAT_concat(:,:,ue,1) )  + Q_2;
        
 end
  if DA==1
      Q_1 =eye(2*N*AP) + kappa^2*(p_transmit(1,ue)/p_pilot(1,ue))*(WdC_0+WdC_1)+ (1-alpha_ADC)/alpha_ADC*((1+kappa^2)*(p_transmit(1,ue)/p_pilot(1,ue))*(WdC_0+WdC_1)+sigma^2*eye(2*N*AP)) +Q_1 ;
      Q_2 = eye(2*N*AP) + kappa^2*(p_transmit(1,ue)/p_pilot(1,ue))*(WdC_0+WdC_1)+ (1-alpha_ADC)/alpha_ADC*((1+kappa^2)*(p_transmit(1,ue)/p_pilot(1,ue))*(WdC_0+WdC_1)+sigma^2*eye(2*N*AP)) + Q_2;
  else
     Q_1 =eye(2*N*AP)+ +Q_1 ;
     Q_2 =eye(2*N*AP)+ +Q_2 ;
  end
 for ue=1:UE   
    mmse_combiner_1(:,ue)= pinv(Q_1)*G_HAT1(:,ue)/Norm_factor_true(1,ue);%%;%sqrt(M(ue,1))  Norm_factor_true(1,ue)
    mmse_combiner_2(:,ue)= pinv(Q_2)*G_HAT2(:,ue)/Norm_factor_true(2,ue);%%;% sqrt(M(ue,2))  Norm_factor_true(2,ue)
    mmse_combiner(:,:,ue)=[(mmse_combiner_1(:,ue) ) (mmse_combiner_2(:,ue) )];
     mrc_combiner_1(:,ue)= G_HAT1(:,ue)/sqrt(trace(R_G_HAT_concat(:,:,ue,1))); % sqrt(trace(R_G_HAT_concat(:,:,ue,1)))   norm(G_HAT1(:,ue))
     mrc_combiner_2(:,ue)= G_HAT2(:,ue)/sqrt(trace(R_G_HAT_concat(:,:,ue,2))); % sqrt(trace(R_G_HAT_concat(:,:,ue,2)))    norm(G_HAT2(:,ue))
     mrc_combiner(:,:,ue)= [mrc_combiner_1(:,ue) mrc_combiner_2(:,ue)];
 end
for ue=1:UE
 G_HAT(:,:,ue)=[G_HAT1(:,ue) G_HAT2(:,ue)];
end

else
   g_hat = zeros(N,AP,UE);
   g_LOS = zeros(N,AP,UE);
   Q = zeros(N,N,AP,UE);
   Q1 = zeros(N,N,AP,UE);
   Q2 =zeros(N*AP,N*AP);
   Q22=zeros(N*AP,N*AP);
%    R_G_HAT_concat=zeros(N*AP,N*AP,UE);
%    C_G_HAT_concat=zeros(N*AP,N*AP,UE);
   y1=zeros(N*AP,UE);
   y=zeros(N,AP,UE);
   C_YY = zeros(N,N,AP,UE);
    C_YY1 = zeros(N*AP,N*AP,UE);

for ap=1:AP
 for ue=1:UE
   for ue1=1:UE
        if rem((ue-ue1),tau_p)==0
   
          y(:,ap,ue)= tau_p*sqrt(p_pilot(ap,ue1))*( (G(:,ap,ue1))) +y(:,ap,ue);
    
          C_YY(:,:,ap,ue) = tau_p*p_pilot(ap,ue1)*R(:,:,ap,ue1) + C_YY(:,:,ap,ue);
        else
           y(:,ap,ue)=y(:,ap,ue); 
           C_YY(:,:,ap,ue) =C_YY(:,:,ap,ue) ;
        end
   end
%       if DA==1
          y_RF(:,ap,ue) = sqrtm(alpha_ADC^2*kappa^2*tau_p*Wd_0(:,:,ap))*(randn(N,1)+1i*randn(N,1))*1/sqrt(2);
          y_ADC(:,ap,ue) = sqrtm(tau_p*alpha_ADC*(1-alpha_ADC)*((1+kappa^2)*Wd_0(:,:,ap)+sigma^2*eye(N)))*(randn(N,1)+1i*randn(N,1))*1/sqrt(2);
          y(:,ap,ue)=sqrt(tau_p)*(randn(N,1)+1i*randn(N,1))*1/sqrt(2) +alpha_ADC*y(:,ap,ue)+y_RF(:,ap,ue)+y_ADC(:,ap,ue);
          C_YY(:,:,ap,ue) =eye(N) +C_YY(:,:,ap,ue)+kappa^2*Wd_0(:,:,ap)+(1-alpha_ADC)/alpha_ADC*((1+kappa^2)*Wd_0(:,:,ap)+eye(N)) ;
%       else
%           y(:,ap,ue)=sqrt(tau_p)*(randn(N,1)+1i*randn(N,1))*1/sqrt(2) +alpha_ADC*y(:,ap,ue);
%           C_YY(:,:,ap,ue) =eye(N) +C_YY(:,:,ap,ue);
%       end
         g_hat(:,ap,ue) = (1/alpha_ADC)*sqrt(p_pilot(ap,ue))*R(:,:,ap,ue)*pinv(C_YY(:,:,ap,ue))*y(:,ap,ue) ;
         R_G_HAT(:,:,ap,ue) = tau_p*p_pilot(ap,ue)*R(:,:,ap,ue)*pinv(C_YY(:,:,ap,ue))*(R(:,:,ap,ue))';
%          g_hat(:,ap,ue) = R_G_HAT(:,:,ap,ue)*(randn(N,1)+1i*randn(N,1))*1/sqrt(2);%
         C_G_HAT(:,:,ap,ue) = R(:,:,ap,ue) -  R_G_HAT(:,:,ap,ue);
%          mrc_combiner(:,ap,ue)=  g_hat(:,ap,ue)/norm( g_hat(:,ap,ue));%norm( g_hat(:,ap,ue));%sqrt(trace(R_G_HAT(:,:,ap,ue)));%
         R_G_concat(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = R(:,:,ap,ue);
         R_G_HAT_concat(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = R_G_HAT(:,:,ap,ue);
         C_G_HAT_concat(((ap-1)*N+1):ap*N,((ap-1)*N+1):ap*N,ue) = C_G_HAT(:,:,ap,ue);
         
 end
end

G_HAT =reshape(g_hat,[N*AP,UE]);

for ue=1:UE
    Q2 = p_transmit(ap,ue)*( G_HAT(:,ue)*G_HAT(:,ue)'+ C_G_HAT_concat(:,:,ue))  + Q2;
   
end
if DA==1
 Q2 =eye(N*AP)+ Q2 + kappa^2*(p_transmit(1,ue)/p_pilot(1,ue))*(WdC_0)+ (1-alpha_ADC)/alpha_ADC*((1+kappa^2)*(p_transmit(1,ue)/p_pilot(1,ue))*(WdC_0)+sigma^2*eye(N*AP)) ;
else
 Q2 = Q2 +eye(N*AP);
end
for ue=1:UE
     mmse_combiner(:,ue) = pinv(Q2)*G_HAT(:,ue)/sqrt(M(ue));% Norm_factor_true(ue) sqrt(M(ue))
     mrc_combiner(:,ue) = G_HAT(:,ue)/sqrt(abs(trace(R_G_HAT_concat(:,:,ue))));%sqrt(abs(trace(R_G_HAT_concat(:,:,ue))));%norm(G_HAT(:,ue))^1;
end
% mrc_combiner=0;
H=0;
end