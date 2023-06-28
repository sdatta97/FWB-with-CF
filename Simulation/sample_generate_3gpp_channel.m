%% Network Setup and Parameter Generation
%
% The tutorial demonstrates how to setup a simple layout with multiple receivers, how to adjust
% parameters manually, generate channel coefficients, and how to calculate parameters from
% the data. The channel model class 'qd_builder' generates correlated values for the LSPs. The
% channel builder then uses those values to create coefficients that have the specific properties
% defined in the builder objects. One important question is therefore: Can the same properties which
% are defined in the builder also be found in the generated coefficients? This is an important test
% to verify, if all components of the channel builder work correctly.

%% Channel model setup and coefficient generation
% We first set up the basic parameters.

% function [channel_coeff, channel_delay, H_fr, c] = sample_generate_3gpp_channel(scenario, cc, BW, num_ue, BS_height, UE_height, Nt, Nr, N, max_dist, min_dist)
function [channel_coeff, channel_delay, H_fr, c] = sample_generate_3gpp_channel(scenario, cc, BW, num_BS, num_ue, BS_height, UE_height, Nt, Nr, N, UE_locations, AP_locations)

% cc = 3.5e9;
% num_ue = 1;
% BS_height = 20;

set(0,'defaultTextFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontSize', 18)                        % Default Font Size
set(0,'defaultAxesFontName','Times')                    % Default Font Type
set(0,'defaultTextFontName','Times')                    % Default Font Type
set(0,'defaultFigurePaperPositionMode','auto')          % Default Plot position
set(0,'DefaultFigurePaperType','<custom>')              % Default Paper Type

s = qd_simulation_parameters;                           % Set up simulation parameters
s.show_progress_bars = 1;                               % Show progress bars
s.center_frequency = cc;                            % Set center frequency
s.samples_per_meter = 1;                                % 1 sample per meter
s.use_absolute_delays = 1;                              % Include delay of the LOS path
s.use_3GPP_baseline = 1;

%% Layout and Channel Generation
% We have one transmitter and 250 receiver positions. Each receiver gets a specific channel.
% However, the receivers LSPs will be correlated. The BS useses a 2-element antenna that transmits a
% linear polarized signal and an left-hand circular polarized signal. This will allow us to verify
% the correct functionality for both polarizations.
 
l = qd_layout(s);% Create new QuaDRiGa layout
l.no_tx = num_BS;
l.no_rx = num_ue;                                            % Set number of MTs
%rx_indices = ones(num_ue,1);
% l.tx_position(1:2,:) = [real(AP_locations), imag(AP_locations)]';
l.tx_position(1:2,:) = AP_locations';
l.tx_position(3,:) = BS_height;
% l.randomize_rx_positions( max_dist , UE_height , UE_height , 0, [], min_dist);
l.rx_position(3,:) = UE_height;
l.rx_position(1:2,:) = UE_locations';
% l.rx_position(1:2,:) = [real(UE_locations), imag(UE_locations)]';

l.set_scenario(scenario);    

if (Nt == 1)
    l.tx_array = qd_arrayant('3gpp-mmw', 1, 0, cc, 1, 5, 0.5, 1, 1, 0, 0);                     % Omni-directional  single BS antenna
elseif (Nt == 4)
    l.tx_array = qd_arrayant('3gpp-mmw', 2, 1, cc, 3, 5, 0.5, 1, 1, 0, 0);
elseif (Nt == 64)
    l.tx_array = qd_arrayant('3gpp-mmw', 4, 8, cc, 3, 5, 0.5, 1, 1, 2, 4);                     % URA muli panel BS antenna
elseif (Nt == 256)
    l.tx_array = qd_arrayant('3gpp-mmw', 4, 8, cc, 3, 5, 0.5, 2, 2, 2, 4);
elseif ((Nt == 32) && (cc <= 6e9))
    l.tx_array = qd_arrayant('3gpp-mmw', 4, 4, cc, 3, 5, 0.5, 1, 1, 2, 4);                      % FR1 BS and NCR
else
    error('Invalid Nt. Only 1, 32, 64 or 256 allowed.');
end

if (Nr == 4)
    l.rx_array = qd_arrayant('3gpp-mmw', 2, 1, cc, 3, 5, 0.5, 1, 1, 0, 0);   % UE
elseif ((Nr == 32) && (cc > 6e9))
    l.rx_array = qd_arrayant('3gpp-mmw', 2, 4, cc, 3, 5, 0.5, 1, 2, 2, 4);   %UE
elseif (Nr == 64)
    l.rx_array = qd_arrayant('3gpp-mmw', 4, 8, cc, 3, 5, 0.5, 1, 1, 2, 4);   %NCR
elseif ((Nr == 32) && (cc <= 6e9))
    l.rx_array = qd_arrayant('3gpp-mmw', 4, 4, cc, 3, 5, 0.5, 1, 1, 2, 4);   % FR1 NCR
elseif ((Nr == 2) && (cc <= 6e9))
    l.rx_array = qd_arrayant('3gpp-mmw', 1, 1, cc, 3, 5, 0.5, 1, 1, 2, 4);   % FR1 UE
else
    error('Invalid Nr. Only 2, 4, 32, 64 allowed.');
end



%%
% We set up the scenario and adjust the parameter range. Then, we generate the channel coefficients.
% In the last step, the arrival angles are obtained from the channel coefficients. This uses only
% the linear polarized transmit singal. 

%beamwidth = calc_beamwidth(l.tx_array);

p = l.init_builder;                                     % Initialize builder
%p.plpar = [];                                           % Remove comment to Disable path-loss

for i = 1:size(p,1)
    p(i,1).scenpar.NumClusters = 15;                             % Reduce paths (for faster processing)
    p(i,1).lsp_xcorr = eye(8);                                   % Disable inter-parameter correlation
end

% p.scenpar.NumClusters = 15;                             % Reduce paths (for faster processing)
% p.lsp_xcorr = eye(8); 

p.gen_parameters;                                       % Generate small-scale-fading parameters
c = p.get_channels;                                     % Generate channel coefficients

% channel_coeff = cell(num_ue,1);
% channel_delay = cell(num_ue,1);
channel_coeff = cell(num_ue,num_BS);
channel_delay = cell(num_ue,num_BS);
% for i = 1:num_ue
%     channel_coeff{i} = c(1,i).coeff;                              % Extract amplitudes and phases
%     channel_delay{i} = c(1,i).delay;
% end
for i = 1:num_ue
    for j = 1:num_BS
        % channel_coeff{i,j} = c(i,j).coeff;                              % Extract amplitudes and phases
        % channel_delay{i,j} = c(i,j).delay;
        channel_coeff{i,j} = c((j-1)*num_ue+i).coeff;                              % Extract amplitudes and phases
        channel_delay{i,j} = c((j-1)*num_ue+i).delay;    
    end
end
%channel_coeff = c.coeff;
%channel_delay = c.delay;
%H_fr = c.fr(BW, (-N/2+1:N/2)/N, 1);                     % N = number of subcarriers
% H_fr = cell(num_ue,1);
% for i = 1:num_ue
%     H_fr{i} = c(1,i).fr(BW, (0:N-1)/N, 1);
% end
H_fr = cell(num_ue,num_BS);
for i = 1:num_ue
    for j = 1:num_BS
        % H_fr{i,j} = c(i,j).fr(BW, (0:N-1)/N, 1);
        H_fr{i,j} = c((j-1)*num_ue+i).fr(BW, (0:N-1)/N, 1);
    end
end
%H_fr = c.fr(BW, (0:N-1)/N, 1);
% set(0,'DefaultFigurePaperSize',[14.5 7.7])              % Adjust paper size for plot
% l.visualize([],[],0);                                   % Plot the layout
% view(-33, 60);                                          % Enable 3D view