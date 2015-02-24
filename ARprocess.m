function x = ARprocess(N)
n = 0:1:N;
rxx = 0;
rho = 0.9;

avgN = 1000;
for i = 1:avgN
    e = sqrt(1-rho^2)*randn(1,N);
    x = filter(1,[1 -rho],e);
    
    rxx = rxx+ xcorr(x, 'coeff');
    
end

rxx = rxx*(1/avgN);
l1 = -(N-1):1:N-1;

figure(1)
stem(l1,rxx), axis([-N, N, -0.05, 1.05]);


sxx = fft(rxx,N);
m = -pi:2*pi/(length(sxx)-1):pi;
figure(2)
plot(m,abs(sxx));

figure(3)
psxx = abs(sxx(1:ceil(length(sxx)/2)));
xsxx = 0.5/length(psxx):0.5/length(psxx):0.5;
plot(xsxx, psxx, 'r');

end


