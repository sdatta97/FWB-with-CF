% function rate_dl = rate_analytical(params, plos2, plos)
function rate_dl = rate_analytical(params, plos2, plos, BETA, ricianFactor)
    Kd_mmw = params.numUE;
    Kd = params.numUE+params.numUE_sub6;
    Kmax = Kd;
    Ntx = params.num_antennas_per_gNB;
    % M = params.numGNB;
    M = params.numGNB_sub6;
    K = params.num_sc_sub6; % no of subcarriers
    scs = params.scs_sub6; %sub carrier spacing
    BW = K*scs;
    a = 1;
    b = 1;
    p_d = 1*Kd;
    % tau_p   = Kd;                  % Uplink training interval. (samples)
    % tau_c   = 200;                 % Coherence Interval. (samples)
    % TAU_FAC = (tau_c-tau_p)/tau_c; % pre-log factor
    TAU_FAC = 1;
    %% initialization of c
    C_v = repmat(sqrt(1./(b*Ntx*sum(BETA,2))),[1 Kd]);
    % E = zeros(M,Kd);
    DS = zeros(Kd,1);
    BU = zeros(Kd,1);
    MUI = zeros(Kd,1);
    % UDI = zeros(Kd,1);
    % TQD = zeros(Kd,1);
    N = ones (Kd,1);
    snr_num = zeros(Kd,1);
    snr_den = zeros(Kd,1);
    rate_dl = zeros(1,Kd);
    for k = 1:Kd_mmw
        DS(k) = DS(k) + a*Ntx*sqrt(p_d)*C_v(:,k)'*BETA(:,k);
        BU(k) = BU(k) + a^2*p_d*Ntx*(C_v(:,k).^2)'*((BETA(:,k)).^2);
        for q = 1:Kd_mmw
             if (q~=k)
                % E(:,q) = sqrt(0.5*(b-a^2)*(C_v(:,q).^2))*(randn(M,1) + 1i*randn(M,1));
                MUI(k) = MUI(k) + plos(q)*a^2*Ntx*p_d*((C_v(:,q)).^2)'*(BETA(:,q).*BETA(:,k));
             end
        end
        for q = 1+Kd_mmw:Kd
            % E(:,q) = sqrt(0.5*(b-a^2)*(C_v(:,q).^2))*(randn(M,1) + 1i*randn(M,1));
            MUI(k) = MUI(k) + a^2*Ntx*p_d*((C_v(:,q)).^2)'*(BETA(:,q).*BETA(:,k));
        end
        % TQD(k) = TQD(k) + abs(sqrt(p_d)*reshape(channel_dl(m,k,:),[Ntx,1]).'*reshape(conj(channel_est_dl(m,q,:)),[Ntx,1])*E(m,q))^2;
        % for l = 1:Ku
        %     UDI(k) = UDI(k) + abs(interUEchannel(k,l)*sqrt(p_u*Theta_v(l)))^2;
        % end
        DS(k) = abs(DS(k))^2;
        snr_num(k) = DS(k);
        % snr_den(k) = MUI(k) + UDI(k) + TQD(k) + N(k);
        snr_den(k) = BU(k) + MUI(k) + N(k);
        rate_dl(k) = BW*TAU_FAC*log2(1+snr_num(k)/snr_den(k));
    end
    for k = 1+Kd_mmw:Kd
        DS(k) = DS(k) + a*Ntx*sqrt(p_d)*C_v(:,k)'*BETA(:,k);
        BU(k) = BU(k) + a^2*Ntx*p_d*((C_v(:,k)).^2)'*(BETA(:,k).^2);
        for q = 1:Kd_mmw
            % E(:,q) = sqrt(0.5*(b-a^2)*(C_v(:,q).^2))*(randn(M,1) + 1i*randn(M,1));
            MUI(k) = MUI(k) + plos(q)*a^2*p_d*Ntx*((C_v(:,q)).^2)'*(BETA(:,q).*BETA(:,k));
        end
        for q = 1+Kd_mmw:Kd
             if (q~=k)
                % E(:,q) = sqrt(0.5*(b-a^2)*(C_v(:,q).^2))*(randn(M,1) + 1i*randn(M,1));
                MUI(k) = MUI(k) + a^2*p_d*Ntx*((C_v(:,q)).^2)'*(BETA(:,q).*BETA(:,k));
             end
        end
        % TQD(k) = TQD(k) + abs(sqrt(p_d)*reshape(channel_dl(m,k,:),[Ntx,1]).'*reshape(conj(channel_est_dl(m,q,:)),[Ntx,1])*E(m,q))^2;
        % for l = 1:Ku
        %     UDI(k) = UDI(k) + abs(interUEchannel(k,l)*sqrt(p_u*Theta_v(l)))^2;
        % end
        DS(k) = abs(DS(k))^2;
        snr_num(k) = DS(k);
        % snr_den(k) = MUI(k) + UDI(k) + TQD(k) + N(k);
        snr_den(k) = BU(k) + MUI(k) + N(k);
        rate_dl(k) = BW*TAU_FAC*log2(1+snr_num(k)/snr_den(k));
    end
end