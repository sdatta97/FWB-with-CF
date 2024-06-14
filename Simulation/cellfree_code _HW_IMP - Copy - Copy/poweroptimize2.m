function [b_mk_opt,b_mk_init,rate_opt,rate0,GEE0,GEE_opt]=poweroptimize2(p_data,UE,AP,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C,combiner_trace,N,B,p_pilot,cvx)
         P_MAX= UE*p_data;
         P_MAX1= P_MAX;% P_MAX;%1;
%          P_TC=0;
         error=1;
         c_mk_n0= abs(c_mk_n0);
         c_mk_n1= abs(c_mk_n1);
         c_mk_d0= abs(c_mk_d0);
         c_mk_d1= abs(c_mk_d1);
         b_mk=ones(AP,UE)*(P_MAX1/(sum(sum(abs(combiner_trace)))))^(0.5);
         b_mk=abs(b_mk);
         b_mk_init=b_mk;
         rate0=SINR(UE,AP,beta,b_mk_init,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C);
         C_AP_DOWN_WCH0=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk_init,B,combiner_trace);
         GEE0=20*rate0/(C_AP_DOWN_WCH0);
         v=1;
         GEE_opt1(1)=10^(-5);
  if cvx==1
    while (error>10^(-4))
        
               [rate_initial,N0,N1,D0,D1]= SINR(UE,AP,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C);
               r=rate_initial;
               [c,pz,x1]=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk,B,combiner_trace);
               t=(c);
               r=(r);
                z=abs(sqrt(r)/t);
               for k=1:UE
                    v0(k) = (N0(k)/D0(k));
                    v1(k) = (N1(k)/D1(k));
                    u0(k)= ((N0(k)+D0(k))/(N0(k)));
                    u1(k)= ((N1(k)+D1(k))/(N1(k)));
                    f_0(k) = ((1+v0(k))/log(2) - beta*sqrt(2)*sqrt(u0(k))/2);
                    f_1(k) = ((1+v1(k))/log(2) - beta*sqrt(2)*sqrt(u1(k))/2);
                    y_0(k) = (sqrt(f_0(k)*N0(k))/(N0(k)+D0(k)));
                    y_1(k) = (sqrt(f_1(k)*N1(k))/(N1(k)+D1(k)));
               end

               for k=1:UE
                   for m=1:AP
                    zx(m,k)=b_mk(m,k)^2*abs(combiner_trace(m,k));
                   end
               end
               J(v)=sum(sum(zx));

             cvx_begin
               variables b_mk(AP,UE) beta_mk(AP,UE) 
               expressions N0(UE) N00(UE) N1(UE) N11(UE) D0(UE) D1(UE) obj(UE) pb(AP,UE) power(AP,UE) a

               for ue=1:UE
                  N00(ue)=0;
                  N11(ue)=0;
                  D0(ue)=0;
                  D1(ue)=0;
                  for ap=1:AP
                      N00(ue)=N00(ue) +(b_mk(ap,ue))*c_mk_n0(ap,ue);
                      N11(ue)=N11(ue) +(b_mk(ap,ue))*c_mk_n1(ap,ue);
                      for ue1=1:UE
                         D0(ue) =D0(ue) + b_mk(ap,ue1)^2*c_mk_d0(ap,ue,ue1);
                         D1(ue) =D1(ue) + b_mk(ap,ue1)^2*c_mk_d1(ap,ue,ue1);
                      end
                 end
%                  N0(ue)=N00(ue)^2;
%                  N1(ue)=N11(ue)^2;
                 D0(ue)=D0(ue)+1;
                 D1(ue)=D1(ue)+1;
               end
        
               for k=1:UE 
                   obj(k)=(2*y_0(k)*sqrt(f_0(k))*N00(k)-y_0(k)^2*(N00(k)^2+D0(k))+ 2*y_1(k)*sqrt(f_1(k))*N11(k)-y_1(k)^2*(N11(k)^2+D1(k))) ;
               end
               for m=1:AP
                   for k=1:UE
                       pb(m,k)=z*sqrt(r)*(b_mk(m,k)^2*abs(pz(m,k))+x1);
                   end
               end
               a=sum(obj)-sum(sum(pb));
               for k=1:UE
                   for m=1:AP
                    power(m,k)=b_mk(m,k)^2*abs(combiner_trace(m,k));
                   end
               end
               
               maximize (a) % -sum(sum(pb))
               subject to
               for k=1:UE
                   for m=1:AP
                    b_mk(m,k)>0;
                   end
               end              
               sum(sum(power))<=P_MAX;
               cvx_end
               
               v=v+1;
               a1(v)=a;
               rate_opt1(v)=SINR(UE,AP,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C);
               C_AP_DOWN_WCH_opt(v)=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk,B,combiner_trace);
               GEE_opt1(v)=20*rate_opt1(v)/(C_AP_DOWN_WCH_opt(v));
               error=abs(GEE_opt1(v)-GEE_opt1(v-1))/(GEE_opt1(v));
%                break;
    end
GEE_opt=GEE_opt1(v);
b_mk_opt=b_mk;
rate_opt=rate_opt1(v);
% rate_opt=1;
% GEE_opt=1;
% lambda=100;%*ones(AP,1);
% error=1;
  else
   
   while error>10^(-4)
      [rate_initial,N0,N1,D0,D1]= SINR(UE,AP,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C);
       c_initial=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk,B,combiner_trace);
        GEE_initial=20*rate_initial/c_initial;
        t1=1;
        r=rate_initial;
        [c,pz,x1,c_p]=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk,B,combiner_trace);
        t=(c);
        r=(r);
        z=abs(sqrt(r)/t);
        for k=1:UE
           v0(k) = (N0(k)/D0(k));
           v1(k) = (N1(k)/D1(k));
           u0(k)= ((N0(k)+D0(k))/(N0(k)));
           u1(k)= ((N1(k)+D1(k))/(N1(k)));
           f_0(k) = ((1+v0(k))/log(2) - beta*sqrt(2)*sqrt(u0(k))/2);
           f_1(k) = ((1+v1(k))/log(2) - beta*sqrt(2)*sqrt(u1(k))/2);
           y_0(k) = (sqrt(f_0(k)*N0(k))/(N0(k)+D0(k)));
           y_1(k) = (sqrt(f_1(k)*N1(k))/(N1(k)+D1(k)));
        end

        for k=1:UE
          for m=1:AP
              zx(m,k)=b_mk(m,k)^2*abs(combiner_trace(m,k));
          end
        end
        J(t1)=sum(sum(zx));
        t2=0;
       startP=0;
       endP=1000;
       lambda=endP;       
       while (endP-startP)/endP >10^(-4)
          rate0=SINR(UE,AP,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C);
          c0=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk,B,combiner_trace);
          GEE1=20*rate0/c0;
          b_mk_den=zeros(AP,UE);
          t2=t2+1;
          for m=1:AP
            for k=1:UE
            
                b_mk_num(m,k)=y_0(k)*sqrt(f_0(k))*c_mk_n0(m,k) +  y_1(k)*sqrt(f_1(k))*c_mk_n1(m,k) - 2*y_0(k)^2*(sqrt(N0(k))-b_mk(m,k)*c_mk_n0(m,k))*c_mk_n0(m,k) - 2*y_0(k)^2*(sqrt(N1(k))-b_mk(m,k)*c_mk_n1(m,k))*c_mk_n1(m,k) ;
                for i=1:UE
                   b_mk_den(m,k)=b_mk_den(m,k) + 2*y_0(i)^2*c_mk_d0(m,k,i) + 2*y_1(i)^2*c_mk_d1(m,k,i) ;
                end
                b_mk_den(m,k)=b_mk_den(m,k)+2*(z^2*c_p+lambda)*sqrt(r)/z*combiner_trace(m,k) + 2*y_0(k)^2*c_mk_n0(m,k)^2 + 2*y_1(k)^2*c_mk_n1(m,k)^2;
                b_mk(m,k) = b_mk_num(m,k)/b_mk_den(m,k);
           end
          end
          power=0;
          for k=1:UE
              for m=1:AP
                  power=b_mk(m,k)^2*combiner_trace(m,k)+power;
              end
          end
          if power<P_MAX
           endP=(endP+startP)/2;

           rate1=SINR(UE,AP,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C);
           c1=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk,B,combiner_trace);
           GEE2(t2)=20*(rate1)/c1;
          else
           startP=(endP+startP)/2;
          end
           lambda= (endP+startP)/0.02;
           
        end
        rate_final=SINR(UE,AP,beta,b_mk,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C);
        c_final=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk,B,combiner_trace);
        GEE_final=20*rate_final/c_final;
        t1=t1+1;
        error=abs(GEE_final-GEE_initial)/GEE_initial;
   end
   b_mk_opt=b_mk;
   rate_opt=SINR(UE,AP,beta,b_mk_opt,c_mk_n0,c_mk_n1,c_mk_d0,c_mk_d1,T_C);
   c_opt=Power_consumption_downlink_wch(N,AP,UE,T_C, p_pilot,1,b_mk,B,combiner_trace);
   GEE_opt=20*rate_opt/c_opt;
  end

end
