function p_out = CF_outage_analytical(params,ue_idx)
    p_out = 0.5;
    N_UE_sub6 = params.N_UE_sub6;
    N_UE_mmW = params.N_UE_mmW;
    T = params.r_min_sub6 (1);
    mu = 3.67;
    d0 = 1;
    r0 = min(abs(params.locationsBS_sub6 - params.UE_locations(ue_idx)));
    R = params.coverageRange_sub6;
    H = hypergeom([1,2],1,-0.1);
    hfunc = @(b,N,r1,r2) ((d0^2)/(b^N*(mu*N+2))*(hypergeom([N,N+(2/mu)],N+(2/mu)+1,-1/(b*(1+r2/d0)^(-mu)))/((1+r2/d0)^(-mu*(N+2/mu)))- hypergeom([N,N+(2/mu)],N+(2/mu)+1,-1/(b*(1+r1/d0)^(-mu)))/((1+r1/d0)^(-mu*(N+2/mu))))...
        - (d0^2)/(b^N*(mu*N+1))*(hypergeom([N,N+(1/mu)],N+(1/mu)+1,-1/(b*(1+r2/d0)^(-mu)))/((1+r2/d0)^(-mu*(N+1/mu)))- hypergeom([N,N+(1/mu)],N+(1/mu)+1,-1/(b*(1+r1/d0)^(-mu)))/((1+r1/d0)^(-mu*(N+1/mu)))));
    term1 = @(b1,n1,N1,b2,N2,r1,r2) (b1^n1*hfunc(b1,n1,r1,r2)*((-1)^(N1-n1)*(factorial(N2+N1-n1-1)/factorial(N2-1))*(b2^-1-b1^-1)^(-N2-(N1-n1)))/(factorial(N1-n1)));
    term2 = @(b1,N1,b2,n2,N2,r1,r2) (b2^n2*hfunc(b2,n2,r1,r2)*((-1)^(N2-n2)*(factorial(N1+N2-n2-1)/factorial(N1-1))*(b1^-1-b2^-1)^(-N1-(N2-n2)))/(factorial(N2-n2)));
    jfunc = @(b1,N1,b2,N2,r1,r2) ((sum(term1(b1,1:1:N1,N1,b2,N2,r1,r2)) + sum(term2(b1,N1,b2,1:1:N2,r1,r2)))/(b1^N1*b2^N2));
    
    fun = @(x) (1./x - 1./(x.^2));
    p_out = 0.5 - (1/pi)*integral(fun,0,T);
end