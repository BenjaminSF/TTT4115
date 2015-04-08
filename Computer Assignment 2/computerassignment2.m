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
    for j = 1:8
        figure(9+i);
        if j == i
            subplot(4,2,j),stem(u(j,:), 'r');
        else
            subplot(4,2,j),stem(u(j,:), 'b');
        end
    end

end

% 2c
h1 = dfilt.df2t(c,1);
h2 = dfilt.df2t(1,c);
hCas = dfilt.cascade(h1,h2);
figure(31);
freqz(hCas);

figure(30);
zplane(c,1);
title('Minimum-phase filter because all zeros are inside the unit circle');

% 2d
y_filt = [];
v = zeros(8,100);
for i = 1:8
    y_filt(i,:) = filter(1,c,y(i,:));
    [v(1,:),v(2,:),v(3,:),v(4,:),v(5,:),v(6,:),v(7,:),v(8,:)] = RFB(y_filt(i,:),k);
    for j = 1:8
        figure(19+i);
        if j == i
            subplot(4,2,j),stem(v(j,:), 'r');
        else
            subplot(4,2,j),stem(v(j,:), 'b');
        end
    end
end
