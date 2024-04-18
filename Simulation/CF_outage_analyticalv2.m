function p_out = CF_outage_analyticalv2(params,ue_idx,lambda_BS,lambda_UE)
    p_out = 0.5;
    N = params.num_antennas_per_gNB;
    N_UE_sub6 = params.N_UE_sub6;
    N_UE_mmW = params.N_UE_mmW;
    B = params.Band;
    rmin = params.r_min_sub6 (1);
    pdTrue = GeneralizedGamma(1.37, 0.98, 1.60);
end