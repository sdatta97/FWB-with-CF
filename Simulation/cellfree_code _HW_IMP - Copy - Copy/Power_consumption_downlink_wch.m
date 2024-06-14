function [C_AP_DOWN_WCH,pz,x1,c_p]=Power_consumption_downlink_wch(N,AP,UE,T_C, P_P_UP,P_DATA,BETA_DOWN_P,B,combiner_trace)
%% Initialization
P_CH_UP           = zeros(UE,1);
P_DATA_DOWN       = zeros(UE,1);
P_C_DOWN          = zeros(UE,1);
P_TC_DOWN         = zeros(UE,1);
P_BH_DOWN         = zeros(UE,1);
P_TOTAL_DOWN_WCH  = zeros(UE,1);
C_AP_DOWN_WCH     = zeros(UE,1);
THETA             = zeros(UE,UE,AP);


c_u    = 0.4*ones(UE,1); % User  power amplifier efficiency
P_C_AP = 0.2*ones(AP,1); % AP Circuit Power (Watt)
P_C_UE = 0.1; % User circuit power (Watt)
V_dd   = 3; % (Volts)
L_min  = 0.5*10^(-6); %(meter)
f_cor  = 10^6; % (Hz)
I_0    = 10*10^(-6); % (Ampere)
C_p    = 1*10^(-12); % (Farad)
P_FIX  = 0.875; % Per AP fixed power consumption (Watt)
P_LO   = 0.1;   % Local Oscillation Power Consumption
M_AP   = 750*10^(9);
P_BT   = 0.25*10^(-9); % Backhaul power (Watt)
BETA_UP_P=ones(UE,1);
TAU_UP  = 2*UE;


% Power consumption for Uplink
for ue=1:UE
    
    % Pilot power consumption
    P_CH_UP(ue) = TAU_UP*P_P_UP(ue)*BETA_UP_P(ue)/(T_C*c_u(ue));
    
    % Data power consumption
    for ue1=1:UE
        for ap=1:AP
            THETA(ue1,ue1,ap) =(combiner_trace(ap,ue1));% transpose(Gmean(:,ue1,ap))*conj(Gmean(:,ue1,ap))+trace(R_G_HAT(:,:,ue1,ap));
        end
    end
    for ap=1:AP
        P_DATA_DOWN(ue) = P_DATA_DOWN(ue) + (T_C-TAU_UP)*P_DATA/(T_C*c_u(ue))*BETA_DOWN_P(ap,ue)^2*THETA(ue,ue,ap);
        pz(ap,ue)=(T_C-TAU_UP)*P_DATA/(T_C*c_u(ue))*THETA(ue,ue,ap);%
    end
    
    % Circuit power consumption
    P_C_AP_T = 0;
   
    for ap=1:AP
        P_C_AP_T = P_C_AP_T + N*P_C_AP(ap); %AP circuit power consumption
%         for n=1:N
%             P_ADC_AP = P_ADC_AP + (3*V_dd^2*L_min*(2*B+f_cor))/(10^(-0.1525*b_AP(n,ap)+4.838));          % AP ADC power consumption
%             P_DAC_AP = P_DAC_AP + (1/2)*V_dd*I_0*(2^(b_AP(n,ap))-1) + b_AP(n,ap)*C_p*(2*B+f_cor)*V_dd^2; % AP DAC power consumption
%         end
    end
%     P_DAC_U(ue) = (1/2)*V_dd*I_0*(2^(b_U(ue))-1) + b_U(ue)*C_p*(2*B+f_cor)*V_dd^2;    % UE DAC power consumption
%     P_ADC_U(ue) = (3*V_dd^2*L_min*(2*B+f_cor))/(10^(-0.1525*b_U(ue)+4.838));          % AP ADC power consumption
            
    P_TC_DOWN(ue) = (AP*P_FIX + P_LO + P_C_AP_T + UE*P_C_UE )/UE ; % Circuit power consumption
    P_CE     = (3*B*AP*N*(TAU_UP+N))/(T_C*M_AP); % Channel estimation power consumption
    P_SP     = (3*B*AP*N*(T_C-TAU_UP))/(T_C*M_AP); % Signal processing power consumption
    
    P_C_DOWN(ue) = (P_TC_DOWN(ue)+P_CE+P_SP);
    
    % Backhaul power consumption
%     P_BH_DOWN(ue) = B*SE_TH(ue)*AP*P_BT;
    
    %Total power consumption
%     P_TOTAL_DOWN_WCH(ue) = P_CH_UP(ue) + P_DATA_DOWN(ue) + P_C_DOWN(ue) + P_BH_DOWN(ue);
    C_AP_DOWN_WCH1(ue)    = P_CH_UP(ue) + P_DATA_DOWN(ue) + P_C_DOWN(ue);
end
 C_AP_DOWN_WCH=sum( C_AP_DOWN_WCH1);
 x1=sum(P_CH_UP)+sum(P_C_DOWN);
 c_p=(T_C-TAU_UP)/(T_C*c_u(ue));
end