x=1:365;
y=[ones(1,8),zeros(1,357)];
Y=fftshift(fft(y));
yy = [zeros(1,10),ones(1,8),zeros(1,10)];
doc bar
ff = -10:17;
figure(1);  
subplot(2,1,1);  
bar(ff,yy);  
title("Ideal Frequency Filter, DC to 8 Cycles");
ylim([-.1,1.1]); 
xlabel("frequency"); 
grid on;
subplot(2,1,2);  bar(x,real(Y));  xlabel("day of year");  ylabel("weight");
title("time-domain equivalent of ideal filter");


y2 = calc_filter(8, 2, 1, 365);
x2=x-floor(365/2);
figure(2);  
subplot(2,1,1);
bar(0:27, y2(1:28));
ylim([-.1, 1.1]);
grid on;
xlabel("frequency");
title("Ideal Filter convolved w/ Gaussian");
Y2 = fft(y2);
subplot(2,1,2);
bar(x2, fftshift(real(Y2)));
xlim([-100,100])
title("time domain equivalent");
