%Script which implements the tasks of problem 2 in the Computer
%Assignment 1 in the course TTT4115 Communication theory

bitrate = [0.75 2 5];
rho = 0.9;

for k = 1:length(bitrate)
    fprintf('----------Bitrate: %g ----------\n', bitrate(k));
    %Problem 2a
    sigma_q = sqrt(2^(-2*bitrate(k))*pi*exp(1)/6); %Noise variance derived in hand notes
    f = @(x) (1-rho^2)./((1-rho*exp(-1j*2*pi*x)).*(1-rho*exp(1j*2*pi*x))); %S_X defined as function of frequency,x
    sqF = @(x) sqrt(f(x));
    lagrange = (sigma_q/(sigma_q^2+1))*integral(sqF,-0.5,0.5); %squareroot of lambda

    freq = -0.5:0.01:0.5;
    G = sqrt(sigma_q^2./(lagrange^2.*f(freq)))-(sigma_q^2)./f(freq); %Making G(f) and H(f) from formulas provided in the compendium
    H = sqrt(lagrange^2*f(freq)/sigma_q^2)-lagrange^2;

    fprintf('Squareroot of lagrange: %.4g, sigma_q: %.4g\n', lagrange, sigma_q);
    Sqy = sigma_q^2 * H;    %Power spectral density of the noise in the Wiener filter
    Sxy = f(freq).*G.*H;     %Power spectral density of the signal in the Wiener filter

    
    
    figure(1);
    textG =sprintf('|G(f)|^2 for bit-rate %g', bitrate(k));
    textH =sprintf('|H(f)|^2 for bit-rate %g', bitrate(k));
    subplot(3,3,3*k-2), plot(freq,G), title(textG); %Plotting the |G(f)|^2 filter
    subplot(3,3,3*k-1), plot(freq,H), title(textH); %Plotting the |H(f)|^2 filter

    
    text = sprintf('Power spectral density of signal(blue) and noise(red) for bit-rate: %g', bitrate(k));
    figure(1);
    subplot(3,3,3*k);
    semilogy(freq, Sqy, 'r'), title(text); %Plotting the power spectral density of the noise in log-plot
    hold on 
    semilogy(freq,Sxy, 'b'); %Plotting the power spectral density of the signal in log-plot (In the same plot as the corresponding noise
    hold off
    SNR = 10*log10( sum(Sxy)/sum(Sqy) ); %Computing signal-to-nois ratio
    fprintf('Theoretical SNR: %3.4g\n', SNR);
    %Problem 2b
    Sx = f(freq);
    FSMG = FrSamp([G(51:end) G(1:50)]); 
    FSMH = FrSamp([H(51:end) H(1:50)]);
    
    FSMG_fft = fft(FSMG,length(FSMG));
    
    FSMH_fft = fft(FSMH,length(FSMH));
    titleFSMG = sprintf('Fourier transformed of FSM of |G(f)|^2 for bit-rate = %g',bitrate(k));
    titleFSMH = sprintf('Fourier transformed of FSM of |H(f)|^2 for bit-rate = %g',bitrate(k)); 

    figure(2);
    subplot(1,length(bitrate),k),plot(freq,abs(FSMG_fft)),title(titleFSMG);
    hold on
    subplot(1,length(bitrate),k),plot(freq,[G(51:end) G(1:50)])
    figure(3);
    subplot(1,length(bitrate),k),stem(freq,abs(FSMH_fft)),title(titleFSMH);
    hold on
    subplot(1,length(bitrate),k),plot(freq,[H(51:end) H(1:50)]);
    FSMG_fft = FSMG_fft';
    FSMH_fft = FSMH_fft';
    
    filterEstimate = [FSMG_fft(52:end) FSMG_fft(1:51)].*[FSMH_fft(52:end) FSMH_fft(1:51)];
    filterTheory = H.*G;
    titleHGestimate = sprintf('Estimated H(f) and G(f) in cascade for bit-rate = %g',bitrate(k));
    titleHGtheory = sprintf('Theoretical H(f) and G(f) in cascade for bit-rate = %g',bitrate(k));
    figure(5+k);
    subplot(1,2,1),plot(freq,abs(filterEstimate)),title(titleHGestimate);
    figure(5+k);
    subplot(1,2,2),plot(freq,abs(filterTheory)),title(titleHGtheory);
    
    error = abs(filterTheory)-abs(filterEstimate);
    figure(10);
    plot(freq,error);
    %Problem 2c
   
    x = ARprocess2(1000000); %Simulated signal
    
    uMid = conv(x,FSMG);
    
    step = sqrt(12* sigma_q^2);
    quantLimit = 8;
    %partition = -7:step:7;
    partition = (step/2):step:quantLimit;
    partition = [-fliplr(partition) partition];
    codebook = [partition(1)-(step/2) partition+(step/2)];
    
    [idx, uQuant] = quantiz(uMid,partition,codebook);
    
    
    figure(15);
    subplot(1,length(bitrate),k),hist(uQuant,codebook);
    %Find entropy
    hQuant = hist(uQuant,codebook);
    entropyQuant = 0; 
    for m = 1:length(hQuant)
        if hQuant(m) ~= 0
            probQuant = hQuant(m)/length(uQuant);
            entropyQuant = entropyQuant - probQuant*log2(probQuant);
        end
    end
    fprintf('Entropy: %g\n', entropyQuant);
    
    u = conv(uQuant, FSMH); %Simulated output from system
    delay = 100; %Measured, is probably floor(length(filter)/2) per filter
    uDelay = u(101:length(u)-100);
    uNoise = x - uDelay;
    uSNR = 10*log10( sum(x.^2)/sum(uNoise.^2) );
    fprintf('Simulated SNR: %g\n', uSNR);
    figure(42);
    textSNR = sprintf('Original signal (blue) and residual noise signal (red), bitrate: %g', bitrate(k));
    subplot(length(bitrate),1,k), plot(x, 'b'), title(textSNR);
    ylim([-5 5]);
    hold on;
    subplot(length(bitrate),1,k), plot(uNoise,'r');
    
end
