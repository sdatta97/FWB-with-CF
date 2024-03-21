function j = jfunc(b1,N1,b2,N2,r1,r2,d0,mu)
    j = 0;
    for n1 = 1:N1
        j = j + term1(b1,n1,N1,b2,N2,r1,r2,d0,mu);
    end
    for n2 = 1:N2
        j = j + term2(b1,N1,b2,n2,r1,r2,d0,mu);
    end
    j = j/(b1^N1*b2^N2);
end