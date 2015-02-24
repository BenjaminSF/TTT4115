%Script which implements the tasks of problem 2 in the Computer
%Assignment 1 in the course TTT4115 Communication theory

bitrate = [0.75 2 5];
rho = 0.9;

for k = 1:length(bitrate)
    %Problem 2a
    sigma_q = sqrt(2^(-2*bitrate(k))*pi*exp(1)/6); %Noise variance derived in hand notes
    f = @(x) (1-rho^2)./((1-rho*exp(-1j*2*pi*x)).*(1-rho*exp(1j*2*pi*x))); %S_X defined as function of frequency,x
    sqF = @(x) sqrt(f(x));
    lagrange = (sigma_q/(sigma_q^2+1))*integral(sqF,-0.5,0.5); %squareroot of lambda

    freq = -0.5:0.01:0.5;
    G = sqrt(sigma_q^2./(lagrange^2.*f(freq)))-(sigma_q^2)./f(freq); %Making G(f) and H(f) from formulas provided in the compendium
    H = sqrt(lagrange^2*f(freq)/sigma_q^2)-lagrange^2;

    fprintf('Squareroot of lagrange: %f, sigma_q: %f, %f bit-rate\n', lagrange, sigma_q, bitrate(k));
    Sqy = sigma_q^2 * H;    %Power spectral density of the noise in the Wiener filter
    Sxy = f(freq).*G.*H;     %Power spectral density of the signal in the Wiener filter

    
    
    figure(1);
    textG =sprintf('|G(f)|^2 for bit-rate %g', bitrate(k));
    textH =sprintf('|H(f)|^2 for bit-rate %g', bitrate(k));
    subplot(3,3,3*k-2), plot(freq,G), title(textG); %Plotting the |G(f)|^2 filter
    subplot(3,3,3*k-1), plot(freq,H), title(textH); %Plotting the |H(f)|^2 filter

    
    text = sprintf('Power spectral density of signal(blue) and noise(red) for bit-rate: %g', bitrate(k));
    subplot(3,3,3*k);
    semilogy(freq, Sqy, 'r'), title(text); %Plotting the power spectral density of the noise in log-plot
    hold on 
    semilogy(freq,Sxy, 'b'); %Plotting the power spectral density of the signal in log-plot (In the same plot as the corresponding noise

    SNR = 10*log10( sum(Sxy)/sum(Sqy) ); %Computing signal-to-nois ratio
    fprintf('SNR: %3.4g, %g bit-rate.\n', SNR, bitrate(k));
    %Problem 2b
    Sx = f(freq);
    FSMG = FrSamp([G(51:end) G(1:50)]); 
    FSMH = FrSamp([H(51:end) H(1:50)]);
    
    FSMG_fft = fft(FSMG,length(FSMG));
    
    FSMH_fft = fft(FSMH,length(FSMH));
    titleFSMG = sprintf('Fourier transformed of FSM of |G(f)|^2 for bit-rate = %g',bitrate(k));
    titleFSMH = sprintf('Fourier transformed of FSM of |H(f)|^2 for bit-rate = %g',bitrate(k)); 

    figure(10);
    subplot(1,length(bitrate),k),stem(abs(FSMG_fft)),title(titleFSMG);
    figure(11);
    subplot(1,length(bitrate),k),stem(abs(FSMH_fft)),title(titleFSMH);
    
end