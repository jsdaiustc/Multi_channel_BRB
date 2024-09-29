clc;clear;close all;
y1_Fig10=load('y1_Fig10.mat');  
y2_Fig10=load('y2_Fig10.mat');  
y3_Fig10=load('y3_Fig10.mat');  
y1=y1_Fig10.y1;  
y2=y2_Fig10.y2;  
y3=y3_Fig10.y3;  

L=length(y1);
Fs=1000;            % Sampling frequency
Y=[y1,y2,y3];
shift=25;
%% HT-FFT
F=(0:L-1)*Fs/L;
y1=Y(:,1);
y1=abs(hilbert(y1));
y1=y1-mean(y1);
F_sample=F(F<=16);
f_y1=abs(fft(y1,L))/sqrt(L);
f_y1=f_y1(F_sample<=16);
f_y1=10*log10(f_y1)+shift;
subplot(5,1,1)
plot(F_sample,f_y1)
axis([0,16,0,40])
pbaspect([5 1 1])
title('HT-FFT')
%% ESPRIT
G=200;
range=[44,76];
y1=Y(:,1);
[omega,f_y1]=ESPRIT_SelectOrder(y1,G,Fs,range);
figure(1)
f_y1=10*log10(abs(f_y1))+shift;
subplot(5,1,2)
stem(omega-60,f_y1)
axis([0,16,0,40])
pbaspect([5 1 1])
title('ESPRIT')
%% SAP-ESPRIT
G = 200;%
range=[0,16];
f_buttord=[40,60];  % The passband of the Butterworth filter
y=SAP(Y);
y=Butterworth(y,f_buttord(1),f_buttord(2));
[omega,f_y]=ESPRIT_SelectOrder(y,G,Fs,range);
f_y=10*log10(abs(f_y))+shift;
subplot(5,1,3)
stem(omega,f_y)
axis([0,16,0,40])
pbaspect([5 1 1])
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
subplot(5,1,4)
stem(omega,gamma)
axis([0,16,0,40])
pbaspect([5 1 1])
title('SBL')
%% Proposed
threshold=0.2;      % The threshold of SVD
Y_bar=preprocess_proposed(Y,Fs,threshold,f_buttord);
omega=0:0.2:20;
[gamma,omega]=real_valued_SBL(Y_bar,Fs,omega);
gamma=10*log10(gamma)+shift;
figure(1)
subplot(5,1,5)
stem(omega,gamma)
axis([0,16,0,40])
pbaspect([5 1 1])
title('Proposed')