%% Plot of the autocorrelation function and corresponding power spectral density for the generated signal with rho = 0.9
N = 200;
l = -N:1:N;
rho = 0.9;
Rxx = rho.^abs(l);

%figure;
%stem(l,Rxx),title('Autocorrelation of Gaussian AR(1)-process'),xlabel('l'),ylabel('Rxx')
w = -pi:0.005*pi:pi;

sigma_x = 1;
Sxx = (1-rho^2)./((1-rho*exp(-1i.*w)).*(1-rho*exp(1i.*w)));
%figure;
v = -3:6/(length(Sxx)-1):3;
%plot(v,abs(Sxx)),title('Power Spectral Density of Gaussian AR(1)-process'),xlabel('w'),ylabel('Sxx');


n = 0:1:N;
rxx = 0;
%sxx2 = 0;
%xTot = 0;
avgN = 1000;
for i = 1:avgN
    e = sqrt(1-rho^2)*randn(1,N);
    x = filter(1,[1 -rho],e);
    %xTot = xTot + x;
    rxx = rxx+ xcorr(x, 'coeff');
    %sxx2 = sxx2 + (1/N)*(fft(x)).^2;
end
%sxx2 = sxx2*(1/avgN);
%rxx2 = xcorr(xTot);
rxx = rxx*(1/avgN);
l1 = -(N-1):1:N-1;

figure(1)
subplot(1,2,1),stem(l1,rxx), axis([-N, N, -0.05, 1.05]);
subplot(1,2,2),stem(l,Rxx), axis([-N, N, -0.05, 1.05]);

figure(3)
plot(l1, rxx, 'r', l, Rxx, 'b');

sxx = fft(rxx,N);
m = -pi:2*pi/(length(sxx)-1):pi;
figure(2)
subplot(1,2,1),plot(m,abs(sxx));
subplot(1,2,2),plot(v,abs(Sxx));

figure(4)
pSxx = abs(Sxx(find(max(abs(Sxx)) == Sxx):end));
psxx = abs(sxx(1:ceil(length(sxx)/2)));
xSxx = 0.5/length(pSxx):0.5/length(pSxx):0.5;
xsxx = 0.5/length(psxx):0.5/length(psxx):0.5;
plot(xsxx, psxx, 'r');
hold on;
plot(xSxx, pSxx, 'b');


