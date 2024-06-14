function [Norm_factor] =normalization (AP,UE,N,R,alpha,p,pol,p1,C_G_HAT_concat,kappa,alpha_ADC,WDC_0,WDC_1)
tau_p=UE;
for ap =1:AP
       for ue=1:UE 
         G_00(:,ap,ue) = sqrt(1-alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)));
       end
 end
if pol==2
for ap =1:AP
       for ue=1:UE 
         G_00(:,ap,ue) = sqrt(1-alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)));
         G_01(:,ap,ue)= sqrt(alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)) ); %+g_LOS(:,ap,ue)
         G_10(:,ap,ue) = sqrt(alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)) );
         G_11(:,ap,ue) = sqrt(1-alpha)*(sqrtm(R(:,:,ap,ue))*1/sqrt(2)*(randn(N,1)+1i*randn(N,1)));
       end
 end
             
  G=zeros(2*N,2,AP,UE);
  for ue=1:UE
    for ap=1:AP
        G(:,:,ap,ue) =[G_00(:,ap,ue) G_01(:,ap,ue);G_10(:,ap,ue) G_11(:,ap,ue)];
                       
    end
  end  

   else
   G=zeros(N,AP,UE);
  for ap=1:AP
     for ue=1:UE
                        
       G(:,ap,ue) = G_00(:,ap,ue)/sqrt(1-alpha);
    end
  end
end
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
         y1(:,ap,ue)= tau_p*sqrt(p(ap,ue))*( (G(:,1,ap,ue1))) +y1(:,ap,ue);
         y2(:,ap,ue) = tau_p*sqrt(p(ap,ue))*( (G(:,2,ap,ue1))) + y2(:,ap,ue);
         C_YY_1(:,:,ap,ue) = tau_p^2*p(ap,ue)*( R_G_1(:,:,ap,ue1))  +C_YY_1(:,:,ap,ue) ;
         C_YY_2(:,:,ap,ue) = tau_p^2*p(ap,ue)*( R_G_2(:,:,ap,ue1))  + C_YY_2(:,:,ap,ue);
         
        else
          y1(:,ap,ue)= y1(:,ap,ue);
         y2(:,ap,ue) =  y2(:,ap,ue);
         C_YY_1(:,:,ap,ue) = C_YY_1(:,:,ap,ue) ;
         C_YY_2(:,:,ap,ue) =C_YY_2(:,:,ap,ue);
        
        end
  
    end
    y1(:,ap,ue)= sqrtm(tau_p*eye(2*N))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2) +y1(:,ap,ue);
    y2(:,ap,ue) = sqrtm(tau_p*eye(2*N))*(randn(2*N,1)+1i*randn(2*N,1))*1/sqrt(2) + y2(:,ap,ue);
    C_YY_1(:,:,ap,ue) =tau_p*eye(2*N) +C_YY_1(:,:,ap,ue) ;
    C_YY_2(:,:,ap,ue) =tau_p*  eye(2*N) + C_YY_2(:,:,ap,ue);
    g_hat_1(:,ap,ue) = tau_p*sqrt(p(ap,ue))*R_G_1(:,:,ap,ue)*pinv(C_YY_1(:,:,ap,ue))*y1(:,ap,ue);
    g_hat_2(:,ap,ue) = tau_p*sqrt(p(ap,ue))*R_G_2(:,:,ap,ue)*pinv(C_YY_2(:,:,ap,ue))*y2(:,ap,ue); 
    g_hat(:,:,ap,ue)=[(g_hat_1(:,ap,ue) ) (g_hat_2(:,ap,ue) )];
%     mrc_combiner(:,:,ap,ue)=[mrc_combiner1(:,ap,ue)  mrc_combiner2(:,ap,ue) ];

  
 end
end
% R_G_HAT(:,:,:,:,1)= R_G_HAT_1;
% R_G_HAT(:,:,:,:,2)= R_G_HAT_2;
G_HAT1=reshape(g_hat_1,[2*N*AP,UE]);
G_HAT2=reshape(g_hat_2,[2*N*AP,UE]);

 for ue=1:UE
      
         Q_1 = p1(1,ue)*( G_HAT1(:,ue)*G_HAT1(:,ue)' +C_G_HAT_concat(:,:,ue,1) +G_HAT2(:,ue)*G_HAT2(:,ue)'+C_G_HAT_concat(:,:,ue,2)) + (kappa^2 + (1-alpha_ADC)*(1+kappa^2)/alpha_ADC)*(WDC_0+WDC_1)  + Q_1;
         Q_2 = p1(1,ue)*( G_HAT2(:,ue)*G_HAT2(:,ue)'+C_G_HAT_concat(:,:,ue,2)+ G_HAT1(:,ue)*G_HAT1(:,ue)' +C_G_HAT_concat(:,:,ue,1)) + (kappa^2 + (1-alpha_ADC)*(1+kappa^2)/alpha_ADC)*(WDC_0+WDC_1)  + Q_2;
        
 end
    Q_1 =eye(2*N*AP)*(1+(1-alpha_ADC)/alpha_ADC) + Q_1 ;
    Q_2 = eye(2*N*AP)*(1+(1-alpha_ADC)/alpha_ADC) + Q_2;
    
 for ue=1:UE   
    Norm_factor_1(ue)=norm(pinv(Q_1)*G_HAT1(:,ue));
    Norm_factor_2(ue)=norm(pinv(Q_2)*G_HAT2(:,ue));
   Norm_factor(:,ue)=[Norm_factor_1(:,ue); Norm_factor_2(:,ue)];
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
   
          y(:,ap,ue)= tau_p*sqrt(p(ap,ue1))*( (G(:,ap,ue1))) +y(:,ap,ue);
    
          C_YY(:,:,ap,ue) = tau_p^2*p(ap,ue1)*R(:,:,ap,ue1) + C_YY(:,:,ap,ue);
        else
           y(:,ap,ue)=y(:,ap,ue); 
           C_YY(:,:,ap,ue) =C_YY(:,:,ap,ue) ;
        end
   end
          y(:,ap,ue)=sqrt(tau_p)*(randn(N,1)+1i*randn(N,1))*1/sqrt(2) +y(:,ap,ue);
          C_YY(:,:,ap,ue) =tau_p*eye(N) +C_YY(:,:,ap,ue) ;
    
         g_hat(:,ap,ue) = tau_p*sqrt(p(ap,ue))*R(:,:,ap,ue)*pinv(C_YY(:,:,ap,ue))*y(:,ap,ue) ;
                        
 end
end

G_HAT =reshape(g_hat,[N*AP,UE]);

for ue=1:UE
    Q2 = p1(ap,ue)*( G_HAT(:,ue)*G_HAT(:,ue)'+ C_G_HAT_concat(:,:,ue))  + Q2;

end 
Q2 =eye(N*AP)+ Q2 ;


for ue=1:UE
    Norm_factor(ue) = norm(pinv(Q2)*G_HAT(:,ue));

end
% mrc_combiner=0;
H=0;
end