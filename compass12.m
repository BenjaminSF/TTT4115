%Kommunikasjonsteori - semesteroppgave 1
%Oppgave 2
bitrate = [0.75 2 5];
for k = 1:length(bitrate)
    rho = 0.9;
    sigma_q = sqrt(2^(-2*bitrate(k))*pi*exp(1)/6);
    f = @(x) (1-rho^2)./((1-rho*exp(-1j*2*pi*x)).*(1-rho*exp(1j*2*pi*x))); %S_X
    sqF = @(x) sqrt(f(x));
    lagrange = (sigma_q/(sigma_q^2+1))*integral(sqF,-0.5,0.5); %Roten av lambda

    frek = -0.5:0.01:0.5;
    G = sqrt(sigma_q^2./(lagrange^2.*f(frek)))-(sigma_q^2)./f(frek);
    H = sqrt(lagrange^2*f(frek)/sigma_q^2)-lagrange^2;

    fprintf('roten av lagrange: %f, sigma_q: %f, %f bitrate\n', lagrange, sigma_q, bitrate(k));
    Sqy = sigma_q^2 * H;
    Sxy = f(frek).*G.*H;
    
    figure(1);
    textG =sprintf('|G(f)|^2 for bitrate %g', bitrate(k));
    textH =sprintf('|H(f)|^2 for bitrate %g', bitrate(k));
    subplot(3,3,3*k-2), plot(frek,G), title(textG);
        %figure(2);
    subplot(3,3,3*k-1), plot(frek,H), title(textH);

%Sqy = sigma_q*(0.998.*f(frek).^(0.5)-0.00139);
%Sxy = 
    %figure(2+k);
    text = sprintf('Bitrate: %g', bitrate(k));
    subplot(3,3,3*k);
    semilogy(frek, Sqy, 'r'), title(text);
    hold on 
    semilogy(frek,Sxy, 'b');

    SNR = 10*log10( sum(Sxy)/sum(Sqy) );
    fprintf('SNR: %3.4g, %g bitrate.\n', SNR, bitrate(k));
    %Oppgave 2B...
    Sx = f(frek);
end