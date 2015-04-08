close all;
clear all;
ener = [1 zeros(1,99)];
nuller = zeros(1,100);
s = zeros(8,100);
k = tan([0.39,0.41,0.45,0.49]*pi);
x = [];
X = [];
orthogonalCheck = [];

%% Problem 1
for i = 1:8
    s(i,:) = ener;
    for j = 1:8
        if j ~= i
            s(j,:) = nuller;
        end
    end
    x(i,:) = TFB(s(1,:),s(2,:),s(3,:),s(4,:),s(5,:),s(6,:),s(7,:),s(8,:),k);
    figure(1);
    subplot(4,2,i),plot(abs(x(i,:)));

    figure(3);
    hold on;
    plot(x(i,:));
end
hold off;

for i = 1:8
    for j = 1:8
        orthogonalCheck(i,j) = dot(x(i,:),x(j,:));
        %figure;
        %plot(ort);
    end
end
fprintf('For dot value equals 0 the channel i and j is orthogonal\n');
orthogonalCheck

for i = 1:8
   X(i,:) = fft(x(i,:));
   figure(5);
   hold on;
   plot(abs(X(i,:)));
end
hold off;

%% Problem 2

c = [0.3336, 0.2975, 0.1328, 0.0729,0.0389,0.0274,0.0172,0.014,0.0098,0.0087,0.0064,0.0061,0.0047,0.0048,0.0037,0.0042,0.0029,0.0046,0.001,0.0086];
% 2a
y = [];
Y = [];
for i = 1:8
    y(i,:) = filter(c,1,x(i,:));
    Y(i,:) =fft(y(i,:));
    figure(6);
    subplot(4,2,i),plot(abs(Y(i,:)));
    figure(7);
    hold on;
    plot(abs(Y(i,:)));

end
hold off;
z = filter(c,1,ener);
Z = fft(z);
figure(8);
subplot(2,1,1),plot(z);
subplot(2,1,2),plot(abs(Z));
% 2b
u = zeros(8,100);

for i = 1:8
    [u(1,:),u(2,:),u(3,:),u(4,:),u(5,:),u(6,:),u(7,:),u(8,:)] = RFB(y(i,:),k);
%     for j = 1:8
%         figure(9+i);
%         if j == i
%             subplot(4,2,j),stem(u(j,:), 'r');
%         else
%             subplot(4,2,j),stem(u(j,:), 'b');
%         end
%     end

end

% 2c
h1 = dfilt.df2t(c,1);
h2 = dfilt.df2t(1,c);
hCas = dfilt.cascade(h1,h2);
freqz(hCas);
checkString = char('is not', 'is');
%fprintf('The cascade %s an all-pass filter\n', checkString(1+ isallpass(hCas),:));

figure(30);
zplane(c,1);
title('Minimum-phase filter because all zeros are inside the unit circle');
fprintf('The filter %s a min-phase filter\n', strtrim(checkString(1+ isminphase(h1),:)));

% 2d
y_filt = [];
v = zeros(8,100);
for i = 1:8
    y_filt(i,:) = filter(1,c,y(i,:));
    [v(1,:),v(2,:),v(3,:),v(4,:),v(5,:),v(6,:),v(7,:),v(8,:)] = RFB(y_filt(i,:),k);
%     for j = 1:8
%         figure(19+i);
%         if j == i
%             subplot(4,2,j),stem(v(j,:), 'r');
%         else
%             subplot(4,2,j),stem(v(j,:), 'b');
%         end
%     end
end

% 2e
wgn = 0.01 .* randn(size(y));
y_wgn = y + wgn;
y_wgn_filt = [];
vn = zeros(8,100);
for i = 1:8
    y_wgn_filt(i,:) = filter(1,c,y_wgn(i,:));
    [vn(1,:),vn(2,:),vn(3,:),vn(4,:),vn(5,:),vn(6,:),vn(7,:),vn(8,:)] = RFB(y_wgn_filt(i,:),k);
%     for j = 1:8
%         figure(39+i);
%         if j == i
%             subplot(4,2,j),stem(vn(j,:), 'r');
%         else
%             subplot(4,2,j),stem(vn(j,:), 'b');
%         end
%     end
end

% 2f
sSin = sin(0.2:0.2:20);
uSin = zeros(8,100);

xSin = TFB(nuller, sSin, nuller, nuller, nuller, nuller, nuller, nuller, k);
ySin = filter(c,1,xSin);
%ySin_filt = filter(1,c,ySin);
[uSin(1,:),uSin(2,:),uSin(3,:),uSin(4,:),uSin(5,:),uSin(6,:),uSin(7,:),uSin(8,:)] = RFB(ySin, k);
figure(50);
for i = 1:8
    if i == 2
        subplot(4,2,i), stem(uSin(i,:), 'r');
    else
        subplot(4,2,i), stem(uSin(i,:), 'b');
    end
end

%% Problem 3
% 3a
theta = 0:0.1:100;
[Yrow, Ycol] = size(Y);
centerFreq = zeros(1,8);
for i = 1:8
    temp = abs(Y(i,1:ceil(Ycol/2)));
    centerFreq(i) = find(temp == max(temp), 1, 'first');
end
centerFreq = pi .* centerFreq ./ Ycol;
chanAmp = zeros(1,8);
[ampH1, wH1] = freqz(h1);

C = zeros(8,length(theta));
L = zeros(8,length(theta));

for i = 1:8
    tmpPos = find(wH1 >= centerFreq(i), 1, 'first');
    chanAmp(i) = (abs(ampH1(tmpPos)))^2;
    C(i,:) = 0.5 * log2(theta * chanAmp(i));  
end
L = round(2.^C);
%figure(301);
for i = 1:8
    figure(302);
    subplot(4,2,i), plot(theta, C(i,:));
    figure(301);
    subplot(4,2,i), plot(theta, L(i,:));
end
% Ltot = sum(L,1);
% figure(303);
% plot(theta, Ltot);

% 3b
S_x = zeros(8,length(theta));
for i = 1:8
    S_x(i,:) = theta - 1/chanAmp(i);
end
figure(311);
for i = 1:8
    subplot(4,2,i), plot(theta, S_x(i,:), 'b', theta, zeros(length(theta)), 'r');
end

