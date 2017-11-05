% channelizer demo 4

% 30 channel analysis channelizer, 15-to-1 down sample and 30 channel
% synthesis filter, 1-to-15 upsample

% analysis filter, Nyquist
h1=sinc(-(4-1/30):1/30:(4-1/30)).*kaiser(239,8)';

h1x=h1.*exp(j*2*pi*(-119:119)*1/30);
h1y=h1.*exp(j*2*pi*(-119:119)*2/30);

% synthesis filter, Pass Nyquist Spectrum, reject first Nyquist Zone at 2. 
h2=remez(238,[0 0.703 1.28 15]/15,{'myfrf',[1 1 0 0]},[1 1.0]);
figure(1)
subplot(2,1,1)
plot(-(4-1/30):1/30:(4-1/30),h1,'linewidth',2)
grid on
axis([-4 4 -0.3 1.2])
title('Impulse Response, 239-Tap Prototype Nyquist Filter for 30-Path Analysis Channelizer, 8-Taps per Path')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(h1/sum(h1),2000)))),'linewidth',2)
hold on
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(h1x/sum(h1),2000)))),'linewidth',2)
plot([+0.4 +0.6],[-1 -1]*6.02,'r','linewidth',2)
plot([0.5 0.5],[-12 -0],'r','linewidth',2)
hold off
grid on
axis([-3 3 -100 10])
title('Frequency Response, Prototype Nyquist Filter for 30-Path Analysis Channelizer')
xlabel('Frequency')
ylabel('Log Magnitude (dB)')

figure(2)
subplot(4,1,1)
plot(-(4-1/30):1/30:(4-1/30),h1,'linewidth',2)
grid on
axis([-4 4 -0.3 1.2])
title('Impulse Response, 239-Tap Prototype Nyquist Filter for 30-Path Analysis Channelizer, 8-Taps per Path')
xlabel('Time Index')
ylabel('Amplitude')

subplot(4,1,2)
plot(-(4-1/30):1/30:(4-1/30),h2/max(h2),'r','linewidth',2)
grid on
axis([-4 4 -0.3 1.2])
title('Impulse Response, 239-Tap Prototype Nyquist Filter for 30-Path Synthesis Channelizer, 8-Taps per Path')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(h1/sum(h1),2000)))),'linewidth',2)
hold on
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(h1y/sum(h1),2000)))),'--','linewidth',2)
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(conj(h1y)/sum(h1),2000)))),'--','linewidth',2)
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(h2,2000)))),'r','linewidth',2)
hold off
grid on
axis([-3 3 -100 10])
title('Frequency Response, Prototype Nyquist Filter for 30-Path Analysis Channelizer (Blue) and for 30-Path Synthesis Channelizer (Red)')
xlabel('Frequency')
ylabel('Log Magnitude (dB)')

% Full bandwidth Impulse Response

x0=[1 zeros(1,500)];
for cc=1:15
x2=zeros(1,500);
hh1=reshape([0 h1],30,8);
hh2=reshape([0 h2],30,8);

m1=1;  % input index offset
m2=1;  % output index offset

v1=zeros(1,30)';
v2=zeros(1,30)';
v3=zeros(1,30)';
reg1=zeros(30,16);

u1=zeros(1,15)';
u2=zeros(1,30)';
u3=zeros(1,30)';
reg2=zeros(30,16);

flg1=0;
flg2=0;

m2=1;

for n=1:15:length(x0)-14
    v1(1:15)=flipud(x0(n:n+14).');
    v1(16:30)=v1(1:15);
    reg1=[v1 reg1(:,1:15)];
    
    for k=1:15
        v2(k)=reg1(k,1:2:16)*hh1(k,:)';
        v2(k+15)=reg1(k+15,2:2:16)*hh1(k+15,:)';
    end
    
    if flg1==0
        flg1=1;
    else
        flg1=0;
        v2=[v2(16:30);v2(1:15)];
    end
    
    v3=30*ifft(v2);
    
    u3=v3;
    u2=15*ifft(u3);
    
    if flg2==0
        flg2=1;
    else
        flg2=0;
        u2=[u2(16:30);u2(1:15)];
    end
    
    reg2=[u2 reg2(:,1:15)];
    
     for k=1:15
        u_tp=reg2(k,1:2:16)*hh2(k,:)';
        u_bt=reg2(k+15,2:2:16)*hh2(k+15,:)';
        u1(k)=u_tp+u_bt;
    end
    x2(m2:m2+14)=u1.';
    m2=m2+15;
end


figure(3)
subplot(3,1,1)
plot(real(x2),'linewidth',2)
grid on
axis([0 450 -0.1 1.1])
title('Impulse Response, Cascade 15-to-1 30-Path Analysis Channelizer and 1-to-15 30-Path Synthesis Channelizer')
ylabel('Amplitude')

subplot(3,1,2)
plot(real(x2),'linewidth',2)
grid on
axis([0 450 -0.0002 .0002])
title('Impulse Response, Zoom to Low level Artifacts')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,3)
plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(x2,1000)))),'linewidth',2)
grid on
title('Spectrum of Full Bandwidth Impulse Response, In-Band Ripple of Perfect Reconstruction Filter Bank')
xlabel('Frequency')
ylabel('Log Mag (dB)')
pause(0.5)
x0=[0 x0];
end

pause

% Fractional Bandwidth Impulse Response
x0=[1 zeros(1,500)];
for cc=1:15
x2=zeros(1,500);
hh1=reshape([0 h1],30,8);
hh2=reshape([0 h2],30,8);

m1=1;  % input index offset
m2=1;  % output index offset

v1=zeros(1,30)';
v2=zeros(1,30)';
v3=zeros(1,30)';
reg1=zeros(30,16);

u1=zeros(1,15)';
u2=zeros(1,30)';
u3=zeros(1,30)';
reg2=zeros(30,16);

flg1=0;
flg2=0;

m2=1;

for n=1:15:length(x0)-14
    v1(1:15)=flipud(x0(n:n+14).');
    v1(16:30)=v1(1:15);
    reg1=[v1 reg1(:,1:15)];
    
    for k=1:15
        v2(k)=reg1(k,1:2:16)*hh1(k,:)';
        v2(k+15)=reg1(k+15,2:2:16)*hh1(k+15,:)';
    end
    
    if flg1==0
        flg1=1;
    else
        flg1=0;
        v2=[v2(16:30);v2(1:15)];
    end
    
    v3=30*ifft(v2);
    
    u3(1:11)=v3(1:11);
    u3(21:30)=v3(21:30);
    u2=15*ifft(u3);
    
    if flg2==0
        flg2=1;
    else
        flg2=0;
        u2=[u2(16:30);u2(1:15)];
    end
    
    reg2=[u2 reg2(:,1:15)];
    
     for k=1:15
        u_tp=reg2(k,1:2:16)*hh2(k,:)';
        u_bt=reg2(k+15,2:2:16)*hh2(k+15,:)';
        u1(k)=u_tp+u_bt;
    end
    x2(m2:m2+14)=u1.';
    m2=m2+15;
end


figure(4)
subplot(3,1,1)
plot(real(x2),'linewidth',2)
grid on
axis([0 450 -0.2  0.8])
title('Impulse Response, Cascade 15-to-1 30-Path Analysis Channelizer and 1-to-15 30-Path Synthesis Channelizer')
ylabel('Amplitude')

subplot(3,1,2)
plot(real(x2),'linewidth',2)
grid on
axis([0 450 -0.001 .001])
title('Impulse Response, Zoom to Low level Artifacts')
xlabel('Time Index')
ylabel('Amplitude')

g=remez(238,[0 10.15 10.85 15]/15,{'myfrf',[1 1 0 0]},[1 1]);
subplot(3,1,3)
plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(x2,1000)))),'linewidth',2)
% hold on
% plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(g,1000)))),'r','linewidth',2)
% hold off
grid on
axis([-15 15 -100 10])
title('Spectrum of Full Bandwidth Impulse Response, In-Band Ripple of Perfect Reconstruction Filter Bank')
xlabel('Frequency')
ylabel('Log Mag (dB)')
pause(0.5)
x0=[0 x0];
end

figure(5)
subplot(4,1,1)
plot(0:499,real(x2(1:500)),'r','linewidth',2)
grid on
axis([0 450 -0.2 0.8])
title('Impulse Response Cascade Polyphase Filter Banks')
ylabel('Amplitude')

subplot(4,1,2)
plot(0:238,g,'linewidth',2)
grid on
axis([0 450 -0.2 0.8])
title('Impulse Response 239-Tap FIR Filter')
ylabel('Amplitude')
xlabel('Time Index')

subplot(4,1,3)
plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(x2,1000)))),'r','linewidth',2)
grid on
axis([-15 15 -100 10])
title('Frequency Response, Cascade Polyphase Filter Bank')
ylabel('Log Mag (dB)')

subplot(4,1,4)
plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(g,1000)))),'linewidth',2)
grid on
axis([-15 15 -100 10])
title('Frequency Response, FIR Filter')
ylabel('Log Mag (dB)')
xlabel('Frequency')
