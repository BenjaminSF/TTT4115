

Ts = 4;

x1 = ones(1,Ts);
x2 = [1 1 -1 -1];
x3 = [1 -1 -1 1];

h1 = ones(1,Ts);
h2 = [-1 -1 1 1];
h3 = [1 -1 -1 1];

titlex1 = sprintf('Plot for x_1(t) when Ts = %g',Ts);
titlex2 = sprintf('Plot for x_2(t) when Ts = %g',Ts);
titlex3 = sprintf('Plot for x_3(t) when Ts = %g',Ts);
titleh1 = sprintf('Plot for h_1(t) when Ts = %g',Ts);
titleh2 = sprintf('Plot for h_2(t) when Ts = %g',Ts);
titleh3 = sprintf('Plot for x_3(t) when Ts = %g',Ts);
titley1 = sprintf('Plot for y_1(t) when Ts = %g',Ts);
titley2 = sprintf('Plot for y_2(t) when Ts = %g',Ts);
titley3 = sprintf('Plot for y_3(t) when Ts = %g',Ts);

figure(1);
subplot(3,3,1),stem(x1),title(titlex1);
subplot(3,3,2),stem(x2),title(titlex2);
subplot(3,3,3),stem(x3),title(titlex3);

figure(1);
subplot(3,3,4),stem(h1),title(titleh1);
subplot(3,3,5),stem(h2),title(titleh2);
subplot(3,3,6),stem(h3),title(titleh3);

y1= conv(h1,x1);
y2 = conv(h2,x2);
y3 = conv(h3,x3);

figure(1);
subplot(3,3,7),stem(y1),axis([0,8,-2,5]),title(titley1);
subplot(3,3,8),stem(y2),axis([0,8,-2,5]),title(titley2);
subplot(3,3,9),stem(y3),axis([0,8,-2,5]),title(titley3);

fprintf('Max value of y1(t) when Ts = %g is: %g\nMax value of y2(t)when Ts = %g is: %g \nMax value of y3(t) when Ts = %g is: %g\n',Ts,max(y1),Ts,max(y2),Ts,max(y3));


