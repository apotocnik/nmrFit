U = rand();
for i=1:2000
    V = rand();
    
    X = 2*sqrt(-2*log(U))*cos(2*pi*V);
    Y = 2*sqrt(-2*log(U))*sin(2*pi*V);
    
    U=V;
    
    D(i) = X;
    F(i) = Y;
end

hist(D,20)
mean(D)
std(D)