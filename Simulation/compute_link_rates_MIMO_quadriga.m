function rate_dl = compute_link_rates_MIMO_quadriga(params,link,ue_idx,sub6ConnectionState)
%% Multi-frequency simulations
%
% This tutorial demonstrates how to perform simultaneous multi-frequency simulations at two carrier
% frequencies: 2.6 GHz and 28 GHz in an Urban-Macrocell deployment. The BS is equipped with two
% different array antennas. A conventional high-gain antenna operates at 2.6 GHz. The higher
% frequency band uses a massive-MIMO array antenna with in an 8x8 dual-polarized setup. The model is
% consistent in both, the spatial domain and the frequency domain. Simulation assumptions are in
% accordance with 3GPP 38.901 v14.1.0 (see Section 7.6.5 Correlation modeling for multi-frequency
% simulations).
%
% Identical parameters for each frequency:
%
% * LOS / NLOS state must be the same
% * BS and MT positions are the same (antenna element positions are different!)
% * Cluster delays and angles for each multi-path component are the same
% * Spatial consistency of the LSPs is identical
%
% Differences:

%
% * Antenna patterns are different for each frequency
% * Path-loss is different for each frequency
% * Path-powers are different for each frequency
% * Delay- and angular spreads are different
% * num_ue-Factor is different
% * XPR of the NLOS components is different

%% Basic setup
% Multiple frequencies are set in the simulation parameters by providing a vector of frequency
% sample points. A new layout is created with one 25 m high BS positions and 100 MT positions. The
% MTs are placed in accordance with the 3GPP assumptions, where 80% of them are situated indoors at
% different floor levels.
s = qd_simulation_parameters;
% s.center_frequency = [2.6e9, 28e9];                     % Assign two frequencies
s.center_frequency = 2.6e9;                     % Assign frequency

l = qd_layout( s );                                     % New QuaDRiGa layout
l.tx_position = [[params.locationsBS; params.locationsBS_sub6], 25*ones(params.numGNB_sub6,1)]';                              % 25 m BS height
ll = qd_layout( s );                                     % New QuaDRiGa layout
ll.tx_position = [[params.locationsBS; params.locationsBS_sub6], 25*ones(params.numGNB_sub6,1)]';   
% l.no_rx = params.numUE + params.numUE_sub6;                                          % 100 MTs
l.no_rx = params.numUE_sub6;                                          % 100 MTs
ll.no_rx = params.numUE;                                          % 100 MTs

% l.randomize_rx_positions( 200, 1.5, 1.5, 0 );           % Assign random user positions
% l.rx_position(1,:) = l.rx_position(1,:) + 220;          % Place users east of the BS
% l.rx_position = [[params.UE_locations; params.UE_locations_sub6], 1.5*ones(params.numUE+params.numUE_sub6,1)]';
l.rx_position = [[params.UE_locations_sub6], 1.5*ones(params.numUE_sub6,1)]';
ll.rx_position = [[params.UE_locations], 1.5*ones(params.numUE,1)]';
% floor = randi(5,1,l.no_rx) + 3;                         % Set random floor levels
% for n = 1:l.no_rx
%     floor( n ) =  randi(  floor( n ) );
% end
% l.rx_position(3,:) = 3*(floor-1) + 1.5;
% 
indoor_rx = l.set_scenario('3GPP_38.901_UMa_LOS',[],[],0.8);    % Set the scenario
indoor_rx_2 = ll.set_scenario('3GPP_38.901_UMa_LOS',[],[],0.8);    % Set the scenario
% l.rx_position(3,~indoor_rx) = 1.5;                      % Set outdoor-users to 1.5 m height

%% Antenna set-up
% Two different antenna configurations are used at the BS. The 2.6 GHz antenna is constructed from 8
% vertically stacked patch elements with +/- 45 degree polarization. The electric downtilt is set to
% 8 degree. The mm-wave antenna uses 64 dual-polarized elements in a 8x8 massive-MIMO array
% configuration. The antennas are assigned to the BS by an array of "qd_arrayant" objects. Rows
% correspond to the frequency, columns to the BS. There is only 1 BS in the layout. The mobile
% terminal uses a vertically polarized omni-directional antenna for both frequencies.

% a_2600_Mhz  = qd_arrayant( '3gpp-3d',  8, 1, s.center_frequency(1), 6, 8 );
% a_28000_MHz = qd_arrayant( '3gpp-3d',  8, 8, s.center_frequency(2), 3 );
a_2600_Mhz  = qd_arrayant( '3gpp-3d',  params.num_antennas_per_gNB/2, params.num_antennas_per_gNB/2, s.center_frequency,6,8);

l.tx_array = a_2600_Mhz;                           % Set 2.6 GHz antenna
ll.tx_array = a_2600_Mhz;                           % Set 2.6 GHz antenna
% l.tx_array(1,1) = a_2600_Mhz;                           % Set 2.6 GHz antenna
% l.tx_array(2,1) = a_28000_MHz;                          % Set 28 Ghz antenna

% l.rx_array = qd_arrayant('omni');                       % Set omni-rx antenna
aa_2600_Mhz = qd_arrayant( '3gpp-3d',  params.N_UE_sub6, 1, s.center_frequency);
aaa_2600_Mhz = qd_arrayant( '3gpp-3d',  params.N_UE_mmW, 1, s.center_frequency);
l.rx_array = aa_2600_Mhz;
ll.rx_array = aaa_2600_Mhz;
%% Coverage preview
% Next, we create a preview of the antenna footprint. We calculate the map for the two frequencies
% including path-loss and antenna patterns. The first plot is for the 2.6 GHz band.

sample_distance = 5;                                    % One pixel every 5 m
% x_min           = -50;                                  % Area to be samples in [m]
% x_max           = 550;
% y_min           = -300;
% y_max           = 300;
x_min = 0;
x_max = 0;
y_min = 0;
y_max = 0;
rx_height       = 1.5;                                  % Mobile terminal height in [m]
tx_power        = 30;                                   % Tx-power in [dBm] per antenna element
i_freq          = 1;                                    % Frequency index for 2.6 GHz

% Calculate the map including path-loss and antenna patterns
[ map, x_coords, y_coords] = l.power_map_w_bl( '3GPP_38.901_UMa_LOS', 'quick',...
    sample_distance, x_min, x_max, y_min, y_max, rx_height, tx_power, i_freq, link, params);
% [ map, x_coords, y_coords] = l.power_map( '3GPP_38.901_UMa_LOS', 'quick',...
%     sample_distance, x_min, x_max, y_min, y_max, rx_height, tx_power, i_freq );
% 
P_db = 10*log10( sum( map{1}, 4 ) );
[ map, x_coords, y_coords] = ll.power_map_w_bl( '3GPP_38.901_UMa_LOS', 'quick',...
    sample_distance, x_min, x_max, y_min, y_max, rx_height, tx_power, i_freq, link, params);
% [ map, x_coords, y_coords] = l.power_map( '3GPP_38.901_UMa_LOS', 'quick',...
%     sample_distance, x_min, x_max, y_min, y_max, rx_height, tx_power, i_freq );
% 
P_db_mmW = 10*log10( sum( map{1}, 4 ) );
%%
% % For the 28 GHz, we get the complex-valued phases for each antenna element in order
% % to calculate a MRT beamformer that points the towards the ground at coordinates x = 200 m and 
% % y = 100 m.
% 
% tx_power        = 10;                                   % Tx-power in [dBm] per antenna element
% i_freq          = 2;                                    % Frequency index for 28 GHz
% 
% % Calculate the map including path-loss and antenna patterns
% [ map, x_coords, y_coords] = l.power_map( '3GPP_38.901_UMa_LOS', 'phase',...
%     sample_distance, x_min, x_max, y_min, y_max, rx_height, tx_power, i_freq );
% % [ map, x_coords, y_coords] = l.power_map( '3GPP_38.901_UMa_LOS', 'detailed',...
% %     sample_distance, x_min, x_max, y_min, y_max, rx_height, tx_power, i_freq );
% % Calculate MRT beamforming weights
% beam_x = find( x_coords >= 200 , 1 );                   % Point the beam at x = 200 and y = 100
% beam_y = find( y_coords >= 100  , 1 );
% w = conj( map{1}( beam_y, beam_x , 1 ,: ) );            % Precoding weights for a MRT beamformer
% w = w ./ sqrt(mean(abs(w(:)).^2));                      % Normalize to unit power
% 
% % Apply the precoding weights to each pixel on the map and calculate the received power
% P_db = map{1} .* w( ones(1,numel(y_coords)), ones(1,numel(x_coords)),:,: );
% P_db = 10*log10( abs( sum( P_db ,4 ) ).^2 );


%% Generate channel coefficients
% Channel coefficients are generated by calling "l.get_channels". The output is an array of QuaDRiGa
% channel objects. The first dimension corresponds to the MTs (100). The second dimension
% corresponds to the number of BSs (1) and the third dimension corresponds to the number of
% frequencies (2).

c = l.get_channels;
cc = ll.get_channels;
num_bs = params.numGNB_sub6;
% num_ue = l.no_rx;
% num_ue_mmW = params.numUE;
num_ue = l.no_rx + ll.no_rx;
num_ue_mmW = ll.no_rx;

% channel_coeff = cell(num_ue,num_bs);
% channel_delay = cell(num_ue,num_bs);
% for i = 1:num_ue
%     for j = 1:num_bs
%         channel_coeff{i,j} = c(i,j).coeff;                              % Extract amplitudes and phases
%         channel_delay{i,j} = c(i,j).delay;
%     end
% end
BW = params.Band; %bandwidth
N = 1; %number of subcarriers

TAU_FAC = (params.tau_c - params.tau_p)/params.tau_c;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(BW) + noiseFigure;
noiseVariance = db2pow(noiseVariancedBm);

N_AP = params.num_antennas_per_gNB;
N_UE_mmW = params.N_UE_mmW;
N_UE_sub6 = params.N_UE_sub6;
Ntx = params.num_antennas_per_gNB;
p_fac = params.p_fac;
p_d = params.rho_tot;
%channel_coeff = c.coeff;
%channel_delay = c.delay;
%H_fr = c.fr(BW, (-N/2+1:N/2)/N, 1);                     % N = number of subcarriers
H_fr = cell(num_ue,num_bs);
BETA = zeros(num_bs,num_ue);
beta_uc = zeros(num_bs,num_ue);
channel_dl_mmW = zeros(num_bs,num_ue_mmW,Ntx,N_UE_mmW);
channel_dl = zeros(num_bs,num_ue - num_ue_mmW,Ntx,N_UE_sub6);
for i = 1:num_ue
    for j = 1:num_bs
        if (i<=num_ue_mmW)
            H_fr{i,j} = cc(i,j).fr(BW, (0:N-1)/N, 1);
        else
            H_fr{i,j} = c(i-num_ue_mmW,j).fr(BW, (0:N-1)/N, 1);
        end
        BETA(j,i) = (mean(abs(H_fr{i,j}),"all"))^2; 
        if (i<=num_ue_mmW)
            channel_dl_mmW(j,i,:,:) = (H_fr{i,j}).';
        else
            channel_dl(j,i-num_ue_mmW,:,:) = (H_fr{i,j}).';
        end
    end
end
channel_est_dl_mmW = channel_dl_mmW;
channel_est_dl = channel_dl;
D = params.D;
%Prepare array to store the number of APs serving a specficic UE
La = zeros(num_ue,1);
%Prepare cell to store the AP indices serving a specficic UE
Serv = cell(num_ue,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ = cell(num_ue,1);
%Construc the above array and cells
for k = 1:num_ue
    servingAPs = find(D(:,k)==1);
    NoservingAPs = find(D(:,k)==0);
    
    Serv{k} = servingAPs;
    NoServ{k} = NoservingAPs;
    
    La(k) = length(servingAPs);
    beta_uc(:,k) = BETA(:,k).*D(:,k);
end

%% initialization of c
eta_eq = zeros(num_bs,num_ue);
if (num_ue_mmW == 0)
    for m = 1:num_bs
        for k = 1:num_ue
            if ismember(m,Serv{k})
                eta_eq(m,k) = 1./(N_AP*N_UE_sub6*sum(BETA(m,:)));
            end
        end
    end
else
    for m = 1:num_bs
        for k = 1:num_ue
            if ismember(m,Serv{k})
                if ((k<=num_ue_mmW) && (sub6ConnectionState(k) == 1))
%                     eta_eq(m,k) = p_fac./(N_AP*(N_UE_mmW*p_fac*beta_uc(m,1:num_ue_mmW)+N_UE_sub6*sum(beta_uc(m,2:num_ue))));
                    eta_eq(m,k) = p_fac./(N_AP*(N_UE_mmW*p_fac*(beta_uc(m,1:num_ue_mmW).*sub6ConnectionState)+N_UE_sub6*sum(beta_uc(m,2:num_ue))));
                elseif (k>num_ue_mmW)
%                     eta_eq(m,k) = 1./(N_AP*(N_UE_mmW*p_fac*beta_uc(m,1:num_ue_mmW)+N_UE_sub6*sum(beta_uc(m,2:num_ue))));
                    eta_eq(m,k) = 1./(N_AP*(N_UE_mmW*p_fac*(beta_uc(m,1:num_ue_mmW).*sub6ConnectionState)+N_UE_sub6*sum(beta_uc(m,2:num_ue))));
                end
            end
        end
    end
end
D_mmW_mmW = zeros(num_ue_mmW,num_ue_mmW,N_UE_mmW,N_UE_mmW);
D_mmW_sub6 = zeros(num_ue_mmW,num_ue-num_ue_mmW,N_UE_mmW,N_UE_sub6);
D_sub6_mmW = zeros(num_ue-num_ue_mmW,num_ue_mmW,N_UE_sub6,N_UE_mmW);
D_sub6_sub6 = zeros(num_ue-num_ue_mmW,num_ue-num_ue_mmW,N_UE_sub6,N_UE_sub6);
for k = 1:num_ue_mmW
    for q = 1:num_ue_mmW
        for m = 1:num_bs
            D_mmW_mmW(k,q,:,:) = reshape(D_mmW_mmW(k,q,:,:),[N_UE_mmW,N_UE_mmW]) + sqrt(eta_eq(m,q))*reshape(channel_dl_mmW(m,k,:,:),[Ntx,N_UE_mmW])'*reshape(conj(channel_est_dl_mmW(m,q,:,:)),[Ntx,N_UE_mmW]);
        end
    end
    for q = 1:num_ue-num_ue_mmW
        for m = 1:num_bs
            D_mmW_sub6(k,q,:,:) = reshape(D_mmW_sub6(k,q,:,:),[N_UE_mmW,N_UE_sub6]) + sqrt(eta_eq(m,q))*reshape(channel_dl_mmW(m,k,:,:),[Ntx,N_UE_mmW])'*reshape(conj(channel_est_dl(m,q,:,:)),[Ntx,N_UE_sub6]);
        end
    end
end
for k = 1:num_ue-num_ue_mmW
    for q = 1:num_ue_mmW
        for m = 1:num_bs
            D_sub6_mmW(k,q,:,:) = reshape(D_sub6_mmW(k,q,:,:),[N_UE_sub6,N_UE_mmW]) + sqrt(eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[Ntx,N_UE_sub6])'*reshape(conj(channel_est_dl_mmW(m,q,:,:)),[Ntx,N_UE_mmW]);
        end
    end
    for q = 1:num_ue-num_ue_mmW
        for m = 1:num_bs
            D_sub6_sub6(k,q,:,:) = reshape(D_sub6_sub6(k,q,:,:),[N_UE_sub6,N_UE_sub6]) + sqrt(eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[Ntx,N_UE_sub6])'*reshape(conj(channel_est_dl(m,q,:,:)),[Ntx,N_UE_sub6]);
        end
    end
end
DS_mmW = zeros(num_ue_mmW,N_UE_mmW);
MSI_mmW = zeros(num_ue_mmW,N_UE_mmW);
MUI_mmW = zeros(num_ue_mmW,N_UE_mmW);
DS_sub6 = zeros(num_ue-num_ue_mmW,N_UE_sub6);
MSI_sub6 = zeros(num_ue-num_ue_mmW,N_UE_sub6);
MUI_sub6 = zeros(num_ue-num_ue_mmW,N_UE_sub6);

% noise_mmW = abs(sqrt(0.5)*(randn(num_ue_mmW,N_UE_mmW) + 1j*randn(num_ue_mmW,N_UE_mmW))).^2;
% noise_sub6 = abs(sqrt(0.5)*(randn(num_ue-num_ue_mmW,N_UE_sub6) + 1j*randn(num_ue-num_ue_mmW,N_UE_sub6))).^2;
noise_mmW = abs(sqrt(0.5*noiseVariance)*(randn(num_ue_mmW,N_UE_mmW) + 1j*randn(num_ue_mmW,N_UE_mmW))).^2;
noise_sub6 = abs(sqrt(0.5*noiseVariance)*(randn(num_ue-num_ue_mmW,N_UE_sub6) + 1j*randn(num_ue-num_ue_mmW,N_UE_sub6))).^2;
snr_num_mmW = zeros(num_ue_mmW,N_UE_mmW);
snr_den_mmW = zeros(num_ue_mmW,N_UE_mmW);
snr_num_sub6 = zeros(num_ue-num_ue_mmW,N_UE_sub6);
snr_den_sub6 = zeros(num_ue-num_ue_mmW,N_UE_sub6);
rate_dl = zeros(num_ue,1);
for k = 1:num_ue_mmW
    if (sub6ConnectionState(k)==1 || k==ue_idx)
        for n = 1:N_UE_mmW
%             DS_mmW(k,n) = p_d*norm(reshape(D_mmW_mmW(k,k,n,:),[1,N_UE_mmW]))^2;
            DS_mmW(k,n) = p_d*(abs(D_mmW_mmW(k,k,n,n)))^2;
            for nn = 1:N_UE_mmW
                if (nn~=n)
%                     MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,k,nn,:),[1,N_UE_mmW]))^2;
                    MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*(abs(D_mmW_mmW(k,k,n,nn)))^2;
                end
            end
            for q = 1:num_ue_mmW
                if (q~=k && sub6ConnectionState(q)==1)
                  MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,q,n,:),[1,N_UE_mmW]))^2;
                end
            end
            for q = 1:num_ue-num_ue_mmW
               MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_sub6(k,q,n,:),[1,N_UE_sub6]))^2;
            end
            snr_num_mmW(k,n) = DS_mmW(k,n);
%             snr_den_mmW(k,n) = MUI_mmW(k,n) + noise_mmW(k,n);
            snr_den_mmW(k,n) = MSI_mmW(k,n) + MUI_mmW(k,n) + noise_mmW(k,n);
            rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
        end
    end
end
for k = 1:num_ue-num_ue_mmW
    for n = 1:N_UE_sub6
%         DS_sub6(k,n) = p_d*norm(reshape(D_sub6_sub6(k,k,n,:),[1,N_UE_sub6]))^2;
        DS_sub6(k,n) = p_d*(abs(D_sub6_sub6(k,k,n,n)))^2;
        for nn = 1:N_UE_sub6
            if (nn~=n)
%                 MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,k,nn,:),[1,N_UE_sub6]))^2;
                MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*(abs(D_sub6_sub6(k,k,n,nn)))^2;
            end
        end
        for q = 1:num_ue_mmW
            if (q~=k && sub6ConnectionState(q)==1)
              MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_mmW(k,q,n,:),[1,N_UE_mmW]))^2;
            end
        end
        for q = 1:num_ue-num_ue_mmW
           MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,q,n,:),[1,N_UE_sub6]))^2;
        end
        snr_num_sub6(k,n) = DS_sub6(k,n);
%         snr_den_sub6(k,n) = MUI_sub6(k,n) + noise_sub6(k,n);
        snr_den_sub6(k,n) = MSI_sub6(k,n) + MUI_sub6(k,n) + noise_sub6(k,n);
        rate_dl(k+num_ue_mmW) = rate_dl(k+num_ue_mmW) + BW*TAU_FAC*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
    end
end
end