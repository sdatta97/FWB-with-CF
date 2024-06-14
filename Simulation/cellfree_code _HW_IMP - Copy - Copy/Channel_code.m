% function Channel_LOS = CH(alpha,theta_dept,ue,ap,pol)
psi = rand(1,4)*2*pi;
ue=10;
ap=16;
pol=4;
channel_LOS = zeros(pol,ap,ue);
alpha=0.1;
theta_dept=2*pi/19;
d_by_lambda=1/20;
for POL=1:pol
  if rem(POL,2)==0  
 for AP =1:ap
     for UE=1:ue
        channel_LOS(POL,AP,UE) =sqrt(1-alpha)* exp(1i*2*pi*(AP-1)*sin(theta_dept)*d_by_lambda); % channelgainovernoise and rician-factor not included.
     end
 end
  else
    for AP =1:ap
     for UE=1:ue
        channel_LOS(POL,AP,UE) =sqrt(alpha)* exp(1i*2*pi*(AP-1)*sin(theta_dept)*d_by_lambda);
     end
    end 
  end
end