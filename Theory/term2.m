function t2 = term2(b1,N1,b2,n2,N2,r1,r2,d0,mu)
    t2 = ((b2^n2)*hfunc(b2,n2,r1,r2,d0,mu)*(((-1)^(N2-n2))*(factorial(N1+N2-n2-1)/factorial(N1-1))*(((1/b1)-(1/b2))^(-N1-(N2-n2)))/(factorial(N2-n2))));
end