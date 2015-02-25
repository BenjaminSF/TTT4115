function x = ARprocess2(N)
rho = 0.9;
e = sqrt(1-rho^2)*randn(1,N);
x = filter(1,[1 -rho],e);
end


