function h = hfunc(b,N,r1,r2,d0,mu)
    try
        h = (((d0^2)/(b^N*(mu*N+2)))*(hypergeom([N,N+(2/mu)],N+(2/mu)+1,-1/(b*(1+r2/d0)^(-mu)))/((1+r2/d0)^(-mu*(N+2/mu)))- hypergeom([N,N+(2/mu)],N+(2/mu)+1,-1/(b*(1+r1/d0)^(-mu)))/((1+r1/d0)^(-mu*(N+2/mu))))...
        - ((d0^2)/(b^N*(mu*N+1)))*(hypergeom([N,N+(1/mu)],N+(1/mu)+1,-1/(b*(1+r2/d0)^(-mu)))/((1+r2/d0)^(-mu*(N+1/mu)))- hypergeom([N,N+(1/mu)],N+(1/mu)+1,-1/(b*(1+r1/d0)^(-mu)))/((1+r1/d0)^(-mu*(N+1/mu)))));
    catch e
        disp(e)
    end
end