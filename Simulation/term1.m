function t1 = term1(b1,n1,N1,b2,N2,r1,r2,d0,mu)
    t1 = ((b1^n1)*hfunc(b1,n1,r1,r2,d0,mu)*(((-1)^(N1-n1))*(factorial(N2+N1-n1-1)/factorial(N2-1))*(((1/b2)-(1/b1))^(-N2-(N1-n1)))/(factorial(N1-n1))));
end