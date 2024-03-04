function p_out = CF_outage_analytical(params,ue_idx,lambda_BS,lambda_UE)
%     p_out = 0.5;
    ps = params.rho_tot;
    N = params.num_antennas_per_gNB;
    N_UE_sub6 = params.N_UE_sub6;
    N_UE_mmW = params.N_UE_mmW;
    T = params.r_min_sub6 (1);
    mu = 3.67;
    d0 = 1;
    r0 = min(sqrt(sum((params.locationsBS_sub6 - params.UE_locations(ue_idx)).^2)));
    R = params.coverageRange_sub6;
    R_inf = 10^6;
%     H = hypergeom([1,2],1,-0.1);
    hfunc = @(b,N,r1,r2) ((d0^2)/(b^N*(mu*N+2))*(hypergeom([N,N+(2/mu)],N+(2/mu)+1,-1/(b*(1+r2/d0)^(-mu)))/((1+r2/d0)^(-mu*(N+2/mu)))- hypergeom([N,N+(2/mu)],N+(2/mu)+1,-1/(b*(1+r1/d0)^(-mu)))/((1+r1/d0)^(-mu*(N+2/mu))))...
        - (d0^2)/(b^N*(mu*N+1))*(hypergeom([N,N+(1/mu)],N+(1/mu)+1,-1/(b*(1+r2/d0)^(-mu)))/((1+r2/d0)^(-mu*(N+1/mu)))- hypergeom([N,N+(1/mu)],N+(1/mu)+1,-1/(b*(1+r1/d0)^(-mu)))/((1+r1/d0)^(-mu*(N+1/mu)))));
    term1 = @(b1,n1,N1,b2,N2,r1,r2) (b1^n1*hfunc(b1,n1,r1,r2)*((-1)^(N1-n1)*(factorial(N2+N1-n1-1)/factorial(N2-1))*((1/b2)-(1/b1))^(-N2-(N1-n1)))/(factorial(N1-n1)));
    term2 = @(b1,N1,b2,n2,N2,r1,r2) (b2^n2*hfunc(b2,n2,r1,r2)*((-1)^(N2-n2)*(factorial(N1+N2-n2-1)/factorial(N1-1))*((1/b1)-(1/b2))^(-N1-(N2-n2)))/(factorial(N2-n2)));
    jfunc = @(b1,N1,b2,N2,r1,r2) ((sum(term1(b1,1:1:N1,N1,b2,N2,r1,r2)) + sum(term2(b1,N1,b2,1:1:N2,r1,r2)))/(b1^N1*b2^N2));
    Ns = 1+pi*lambda_BS*(R^2-r0^2);
    N_bar = ceil(lambda_UE*pi*R^2);
    F1 = @(t) (exp(-2*pi*lambda_BS*(0.5*(R^2-r0^2) - jfunc(-1j*t*ps,N_bar,1j*(t/T)*ps*Ns,N,r0,R))));
    F2 = @(t) (exp(-2*pi*lambda_BS*(0.5*(R_inf^2-R^2) - jfunc(-1j*t*ps,N_bar,0,N-N_bar,R,R_inf))));
    fun = @(t) ((1/t)*fourier(exp(1j*t*(1-1j*t*ps*(1+r0/d0)^(-mu))^(-N_bar)*(1+1j*(t/T)*ps*Ns*(1+r0/d0)^(-mu))^(-N))*F1(t)*F2(t)));
    syms t
    p_out = 0.5 - (1/pi)*integral(fun(t),0,inf);
end