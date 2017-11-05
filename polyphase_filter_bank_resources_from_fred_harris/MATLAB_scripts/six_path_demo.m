% six path filter, 80 dB dynamic range
% 6-dB BW = +/-10 kHz, < channel spacing 20 kHz, fs = 120 kHz 
% 100 dB BW = +/-15 kHz, df=10 kHZ
N=1.4*(120/10)*(80/22); % N = 61... select 65 

phi=2*pi*(-32:32)*10/120;
h0=sin(phi)./phi;
h0(33)=1;
h0=h0.*kaiser(65,8)';
%h0=round((2^15)*(h0))/(2^15);
%h0=6*remez(64,[0 5 10 10 15 60]/60,[1 1 0.500 0.500 0 0],[1 1 1]);
%h0=6*remez(64,[0 5 15 60]/60,[1 1 0 0],[1 1]);

figure(1)
subplot(3,1,1)
plot(0:64,h0,'linewidth',2)
grid on
axis([-5 70 -0.3 1.2])
set(gca,'fontsize',12)
title('Impulse Respone, 65-Tap Protoype Low Pass FIR Nyqusit Filter for 6-Path Analysis Channelizer')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,2)
f_h0=fftshift(20*log10(abs(fft(h0*20/120,2000))));
plot((-0.5:1/2000:0.5-1/2000)*120,f_h0,'linewidth',2)
hold on
plot([-60 -15 -15],[-80 -80 -20],'r','linewidth',2)
plot([+60 +15 +15],[-80 -80 -20],'r','linewidth',2)
plot([-10 -10 +10 +10],[-80 0 0 -80],'r','linewidth',2)
plot([-10 -10],[-90 10],'--r','linewidth',2)
plot([+10 +10],[-90 10],'--r','linewidth',2)
hold off
grid on
axis([-60 60 -100 10])
set(gca,'fontsize',12)
title('Frequency Response and Design Parameters Spectral Mask')
xlabel('Frequency (MHz)')
ylabel('Log Mag (dB)')

subplot(3,2,5)
plot((-0.5:1/2000:0.5-1/2000)*120,f_h0,'linewidth',2)
hold on
plot([-5 -5 +5 +5],[-0.002 -0.001 -0.001 -0.002],'r','linewidth',2)
plot([-5 -5 +5 +5],[+0.002 +0.001 +0.001 +0.002],'r','linewidth',2)
hold off
grid on
axis([-10 10 -0.002 0.002])
set(gca,'fontsize',12)
title('Zoom to Pass Band Ripple and Design Parameters Spectral Mask')
xlabel('Frequency (MHz)')
ylabel('Log Mag (dB)')

subplot(3,2,6)
plot((-0.5:1/2000:0.5-1/2000)*120,f_h0,'linewidth',2)
hold on
plot([0 +5 +5],[-0.1 -0.1 -0.2],'r','linewidth',2)
plot([0 +5 +5],[+0.1 +0.1 +0.2],'r','linewidth',2)
plot([+20 +15 +15],[-80 -80 -20],'r','linewidth',2)
plot([0 +5 +5],[0 0 -90],'r','linewidth',2)
plot([+10 +10],[-90 5],'--r','linewidth',2)
hold off
grid on
axis([0 20 -100 5])
set(gca,'fontsize',12)
title('Transition Bandwidth Detail and Spectral Mask')
xlabel('Frequency (MHz)')
ylabel('Log Mag (dB)')

aa=0;
g0=remez(64,[0 15-aa 25+aa 60]/60,{'myfrf',[1 1 0 0]},[2 1]);
%g0=remez(64,[0 15-aa 25+aa 60]/60,{'myfrf',[1 1 0 0]},[1 1]);
%g0=round((2^17)*(g0))/(2^17);

figure(2)
subplot(3,1,1)
plot(0:64,g0/max(g0),'linewidth',2)
grid on
axis([-5 70 -0.3 1.2])
set(gca,'fontsize',12)
title('Impulse Respone, 65-Tap Protoype Low Pass Remez FIR Filter for 6-Path Synthesis Channelizer')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,2)
f_g0=fftshift(20*log10(abs(fft(g0,2000))));
plot((-0.5:1/2000:0.5-1/2000)*120,f_g0,'linewidth',2)
hold on
plot([-60 -25 -25],[-80 -80 -20],'r','linewidth',2)
plot([+60 +25 +25],[-80 -80 -20],'r','linewidth',2)
plot([-15 -15 +15 +15],[-80 0 0 -80],'r','linewidth',2)
plot([-20 -20],[-90 10],'--r','linewidth',2)
plot([+20 +20],[-90 10],'--r','linewidth',2)
hold off
grid on
axis([-60 60 -100 10])
set(gca,'fontsize',12)
title('Frequency Response and Design Parameters Spectral Mask')
xlabel('Frequency (MHz)')
ylabel('Log Mag (dB)')

subplot(3,2,5)
plot((-0.5:1/2000:0.5-1/2000)*120,f_g0,'linewidth',2)
hold on
plot([-15 -15 +15 +15],[-0.002 -0.001 -0.001 -0.002],'r','linewidth',2)
plot([-15 -15 +15 +15],[+0.002 +0.001 +0.001 +0.002],'r','linewidth',2)
hold off
grid on
axis([-20 20 -0.002 0.002])
set(gca,'fontsize',12)
title('Zoom to Pass Band Ripple and Design Parameters Spectral Mask')
xlabel('Frequency (MHz)')
ylabel('Log Mag (dB)')

subplot(3,2,6)
plot((-0.5:1/2000:0.5-1/2000)*120,f_g0,'linewidth',2)
hold on
plot([10 +15 +15],[-0.1 -0.1 -0.2],'r','linewidth',2)
plot([10 +15 +15],[+0.1 +0.1 +0.2],'r','linewidth',2)
plot([+30 +25 +25],[-80 -80 -20],'r','linewidth',2)
plot([10 +15 +15],[0 0 -90],'r','linewidth',2)
plot([+20 +20],[-90 5],'--r','linewidth',2)
hold off
grid on
axis([10 30 -100 5])
set(gca,'fontsize',12)
title('Transition Bandwidth Detail and Spectral Mask')
xlabel('Frequency (MHz)')
ylabel('Log Mag (dB)')

%%%%%%%%%%%%%%%%%%%%% 3-to-1 Down sample
hh0=reshape([0 h0],6,11);
gg0=reshape([0 g0],6,11);

reg1=zeros(6,22);
reg2=zeros(6,22);
v0=zeros(1,6)';
v1=zeros(1,6)';
v2=zeros(1,6)';
u3=zeros(1,6)';
flg1=0;
flg2=0;

m1=1;
m2=0;

x0=zeros(1,200);
x0(3)=1;
x1=zeros(6,33);
x2=zeros(1,200);

for n=1:3:200-3
    v0(1:3)=fliplr(x0(n:n+2)).';
    v0(4:6)=v0(1:3);
    
    reg1=[v0 reg1(:,1:21)];
    for k=1:3
        v1(k)=reg1(k,1:2:22)*hh0(k,:)';
        v1(k+3)=reg1(k+3,2:2:22)*hh0(k+3,:)';
    end
    
    if flg1==0
        flg1=1;
    else
        flg1=0;
        v1=[v1(4:6);v1(1:3)];
    end
    
    v2=6*ifft(v1);
    x1(:,m1)=v2;
    m1=m1+1;   %finished analysis filter, 6-Impulse responses in x1
    
    u3=3*ifft(v2);  % start synthesis filter
    if flg2==0
        flg2=1;
    else
        flg2=0;
        u3=[u3(4:6);u3(1:3)];
    end
    
    reg2=[u3 reg2(:,1:21)];
    
    for k=1:3
        p1=reg2(k,1:2:22)*gg0(k,:)';
        p2=reg2(k+3,2:2:22)*gg0(k+3,:)';
        x2(m2+k)=(p1+p2);
    end
    m2=m2+3;
end

figure(3)

for k=1:6
subplot(2,3,k)
plot(0:22,real(x1(k,1:23)),'linewidth',2)
hold on
plot(0:22,imag(x1(k,1:23)),'r','linewidth',2)
hold off
grid on
axis([0 22 -0.3 1.2])
set(gca,'fontsize',12)
title(['Impulse Response, Chan(',num2str(k-1),')'],'fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)
end


figure(4)

for k=1:6
subplot(2,3,k)
plot((-0.5:1/128:0.5-1/128)*40,fftshift(20*log10(abs(fft(x1(k,1:23)/2,128)))),'linewidth',2)
grid on
axis([-20 22 -100 10])
set(gca,'fontsize',12)
title(['Frequency Response, Chan(',num2str(k-1),')'],'fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Magnitude (dB)','fontsize',14)
end

figure(6)
subplot(3,1,1)
plot(0:130,x2(1:131),'linewidth',2)
grid on
axis([0 130 -0.1 1.2])
set(gca,'fontsize',12)
title('Impulse Response, Cascade 6-Path Analysis and 6-path Synthesis Channelizers','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

subplot(3,1,2)
plot(0:130,x2(1:131),'linewidth',2)
grid on
axis([0 130 -0.00001 0.00002])
set(gca,'fontsize',12)
title('Zoom To Impulse Response Low Level Artifacts','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

subplot(3,1,3)
plot((-0.5:1/2048:0.5-1/2048)*120,fftshift(20*log10(abs(fft(x2,2048)))),'linewidth',2)
hold on
plot((-0.5:1/512:0.5-1/512)*40,fftshift(20*log10(abs(fft(x1(k,1:23)/2,512)))),'r','linewidth',2)
plot((-0.5:1/512:0.5-1/512)*40+20,fftshift(20*log10(abs(fft(x1(k,1:23)/2,512)))),'r','linewidth',2)
plot((-0.5:1/512:0.5-1/512)*40-20,fftshift(20*log10(abs(fft(x1(k,1:23)/2,512)))),'r','linewidth',2)
hold off
grid on
axis([-60 60 -0.001 0.001])
set(gca,'fontsize',12)
title('Frequency Response, Zoom to Pass Band Ripple Levels','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Magnitude (dB)','fontsize',14)
