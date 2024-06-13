clc;clear;close all;
y1=load('y1.mat');  
y2=load('y2.mat');  
y3=load('y3.mat');  
y1=y1.y1;  
y2=y2.y2;  
y3=y3.y3;  

L=length(y1);
Fs=1000;            % Sampling frequency
Y=[y1,y2,y3];
shift=25;
%% FFT
y1=Y(:,1);
F=(0:L-1)*Fs/L;
F_sample=F(F<=77&F>=43);
f_y1=abs(fft(y1,L))/sqrt(L);
f_y1=f_y1(F<=77&F>=43);
f_y1=10*log10(f_y1)+shift;
figure(1)
subplot(3,2,1)
plot(F_sample-60,f_y1)
axis([0,16,0,40])
title('FFT')
%% HT-FFT
F=(0:L-1)*Fs/L;
y1=Y(:,1);
y1=abs(hilbert(y1));
y1=y1-mean(y1);
F_sample=F(F<=16);
f_y1=abs(fft(y1,L))/sqrt(L);
f_y1=f_y1(F_sample<=16);
f_y1=10*log10(f_y1)+shift;
subplot(3,2,2)
plot(F_sample,f_y1)
axis([0,16,0,40])
title('HT-FFT')
%% ESPRIT
G=200;
range=[44,76];
y1=Y(:,1);
[omega,f_y1]=ESPRIT_SelectOrder(y1,G,Fs,range);
figure(1)
f_y1=10*log10(abs(f_y1))+shift;
subplot(3,2,3)
stem(omega-60,f_y1)
axis([0,16,0,40])
title('ESPRIT')
%% SAP-ESPRIT
G = 200;%
range=[0,16];
f_buttord=[40,60];  % The passband of the Butterworth filter
y=SAP(Y);
y=Butterworth(y,f_buttord(1),f_buttord(2));
[omega,f_y]=ESPRIT_SelectOrder(y,G,Fs,range);
f_y=10*log10(abs(f_y))+shift;
subplot(3,2,4)
stem(omega,f_y)
axis([0,16,0,40])
title('SAP-ESPRIT')
%% SBL
ia=Y(:,1);
baseline = cos( 1/Fs*2*pi*(1:L)*60 );
ia = ia .* baseline.';
ia = hilbert(ia) - mean( hilbert(ia) );
ia = ia ./ max( abs(ia) );
r = 0.2;
omega = 60 + (40:r:80);
N_all = (0:1:length(ia)-1)';
[mu1,gamma,omega] = VBI_offgrid_CGDP(ia,N_all,omega,Fs);
omega=omega-120;
gamma=10*log10(gamma)+shift;
subplot(3,2,5)
stem(omega,gamma)
axis([0,16,0,40])
title('SBL')
%% Our method
threshold=0.2;      % The threshold of SVD
Y_bar=preprocess_proposed(Y,Fs,threshold,f_buttord);
omega=0:0.2:20;
[gamma,omega]=real_valued_SBL(Y_bar,Fs,omega);
gamma=10*log10(gamma)+shift;
figure(1)
subplot(3,2,6)
stem(omega,gamma)
axis([0,16,0,40])
title('Our method')