% channelizer demo 5a

% 4 channel analysis channelizer, 2-to-1 down sample and 30 channel
% synthesis filter, 1-to-15 upsample

% analysis filter, Nyquist
h1=sinc(-(4-1/4):1/4:(4-1/4)).*kaiser(31,8)';

h1x=h1.*exp(j*2*pi*(-15:15)*1/4);
h1y=h1.*exp(j*2*pi*(-15:15)*2/4);

% sythesis filter, Pass Nyquist Spectrum, reject first Nyquist Zone at 2. 
h2=remez(298,[0 0.703 1.28 15]/15,{'myfrf',[1 1 0 0]},[1 1.0]);
figure(1)
subplot(2,1,1)
plot(-(4-1/4):1/4:(4-1/4),h1,'linewidth',2)
grid on
axis([-4 4 -0.3 1.2])
title('Impulse Response, 31-Tap Prototype Nyquist Filter for 4-Path Analysis Channelizer, 8-Taps per Path')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/2000:0.5-1/2000)*4,fftshift(20*log10(abs(fft(h1/sum(h1),2000)))),'linewidth',2)
hold on
plot((-0.5:1/2000:0.5-1/2000)*4,fftshift(20*log10(abs(fft(h1x/sum(h1),2000)))),'linewidth',2)
plot([+0.4 +0.6],[-1 -1]*6.02,'r','linewidth',2)
plot([0.5 0.5],[-12 -0],'r','linewidth',2)
hold off
grid on
axis([-2 2 -100 10])
title('Frequency Response, Prototype Nyquist Filter for 4-Path Analysis Channelizer')
xlabel('Frequency')
ylabel('Log Magnitude (dB)')

figure(2)
subplot(4,1,1)
plot(-(4-1/4):1/4:(4-1/4),h1,'linewidth',2)
grid on
axis([-4 4 -0.3 1.2])
title('Impulse Response, 31-Tap Prototype Nyquist Filter for 4-Path Analysis Channelizer, 8-Taps per Path')
xlabel('Time Index')
ylabel('Amplitude')

subplot(4,1,2)
plot(-(5-1/30):1/30:(5-1/30),h2/max(h2),'r','linewidth',2)
grid on
axis([-4 4 -0.3 1.2])
title('Impulse Response, 299-Tap Prototype Nyquist Filter for 30-Path Synthesis Channelizer, 10-Taps per Path')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/2000:0.5-1/2000)*4,fftshift(20*log10(abs(fft(h1/sum(h1),2000)))),'linewidth',2)
hold on
plot((-0.5:1/2000:0.5-1/2000)*4,fftshift(20*log10(abs(fft(h1y/sum(h1),2000)))),'--','linewidth',2)
plot((-0.5:1/2000:0.5-1/2000)*4,fftshift(20*log10(abs(fft(conj(h1y)/sum(h1),2000)))),'--','linewidth',2)
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(h2,2000)))),'r','linewidth',2)
hold off
grid on
axis([-2 2 -100 10])
title('Frequency Response, Prototype Nyquist Filter for 4-Path Analysis Channelizer (Blue) and for 30-Path Synthesis Channelizer (Red)')
xlabel('Frequency')
ylabel('Log Magnitude (dB)')

% Full bandwidth Impulse Response
x0=[1 zeros(1,520)];
x2=zeros(1,520);
hh1=reshape([0 h1],4,8);
hh2=reshape([0 h2],30,10);

m1=1;  % input index offset
m2=1;  % output index offset

v1=zeros(1,4)';
v2=zeros(1,4)';
v3=zeros(1,4)';
reg1=zeros(4,16);

u1=zeros(1,15)';
u2=zeros(1,30)';
u3=zeros(1,30)';
reg2=zeros(30,20);

flg1=0;
flg2=0;

m2=1;

for n=1:2:length(x0)-1
    v1(1:2)=flipud(x0(n:n+1).');
    v1(3:4)=v1(1:2);
    reg1=[v1 reg1(:,1:15)];
    
    for k=1:2
        v2(k)=reg1(k,1:2:16)*hh1(k,:)';
        v2(k+2)=reg1(k+2,2:2:16)*hh1(k+2,:)';
    end
    
    if flg1==0
        flg1=1;
    else
        flg1=0;
        v2=[v2(3:4);v2(1:2)];
    end
    
    v3=30*ifft(v2);
    
    u3(1:2)=v3(1:2);
    u3(30)=v3(4);
    u2=2*ifft(u3);
    
    if flg2==0
        flg2=1;
    else
        flg2=0;
        u2=[u2(16:30);u2(1:15)];
    end
    
    reg2=[u2 reg2(:,1:19)];
    
     for k=1:15
        u_tp=reg2(k,1:2:20)*hh2(k,:)';
        u_bt=reg2(k+15,2:2:20)*hh2(k+15,:)';
        u1(k)=u_tp+u_bt;
    end
    x2(m2:m2+14)=u1.';
    m2=m2+15;
end


figure(3)
subplot(3,1,1)
plot(10*real(x2),'linewidth',2)
hold on
hold off
grid on
axis([0 520 -0.3 1.1])
title('Impulse Response, Cascade 2-to-1 4-Path Analysis Channelizer and 1-to-15 30-Path Synthesis Channelizer')
ylabel('Amplitude')

subplot(3,1,2)
plot(10*real(x2),'linewidth',2)
grid on
axis([200 325 -0.3 1.1])
title('Impulse Response, Zoom to Main Lobe Response')
xlabel('Time Index')
ylabel('Amplitude')

subplot(3,1,3)
plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(x2,1000)))),'linewidth',2)
grid on
axis([-15 15 -120 10])
title('Spectrum of Full Bandwidth Impulse Response, In-Band Ripple of Perfect Reconstruction Filter Bank')
xlabel('Frequency')
ylabel('Log Mag (dB)')



% Fractional Bandwidth Impulse Response
x0=[1 zeros(1,520)];
x2=zeros(1,520);
hh1=reshape([0 h1],4,8);
hh2=reshape([0 h2],30,10);

m1=1;  % input index offset
m2=1;  % output index offset

v1=zeros(1,4)';
v2=zeros(1,4)';
v3=zeros(1,4)';
reg1=zeros(4,16);

u1=zeros(1,15)';
u2=zeros(1,30)';
u3=zeros(1,30)';
reg2=zeros(30,20);

flg1=0;
flg2=0;

m2=1;

for n=1:2:length(x0)-1
    v1(1:2)=flipud(x0(n:n+1).');
    v1(3:4)=v1(1:2);
    reg1=[v1 reg1(:,1:15)];
    
    for k=1:2
        v2(k)=reg1(k,1:2:16)*hh1(k,:)';
        v2(k+2)=reg1(k+2,2:2:16)*hh1(k+2,:)';
    end
    
    if flg1==0
        flg1=1;
    else
        flg1=0;
        v2=[v2(3:4);v2(1:2)];
    end
    
    v3=30*ifft(v2);
    
    u3(7:9)=[v3(4);v3(1:2)];
    
    u2=2*ifft(u3);
    
    if flg2==0
        flg2=1;
    else
        flg2=0;
        u2=[u2(16:30);u2(1:15)];
    end
    
    reg2=[u2 reg2(:,1:19)];
    
     for k=1:15
        u_tp=reg2(k,1:2:20)*hh2(k,:)';
        u_bt=reg2(k+15,2:2:20)*hh2(k+15,:)';
        u1(k)=u_tp+u_bt;
    end
    x2(m2:m2+14)=u1.';
    m2=m2+15;
end


figure(4)
subplot(3,1,1)
plot(10*real(x2),'linewidth',2)
hold on
plot(10*imag(x2),'r','linewidth',2)
hold off
grid on
axis([0 520 -1.0  1.0])
title('Impulse Response, Cascade 2-to-1 4-Path Analysis Channelizer and 1-to-15 30-Path Synthesis Channelizer')
ylabel('Amplitude')

subplot(3,1,2)
plot(10*real(x2),'linewidth',2)
hold on
plot(10*imag(x2),'r','linewidth',2)
hold off
grid on
axis([200 325 -1.0 1.0])
title('Impulse Response, Zoom to Main Lobe Response')
xlabel('Time Index')
ylabel('Amplitude')

g=remez(298,[0 1.3 1.8 15]/15,{'myfrf',[1 1 0 0]},[1 1]);
gg=g.*exp(j*2*pi*(-149:149)*7./30);
% 
subplot(3,1,3)
plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(x2,1000)))),'linewidth',2)
% hold on
% plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(gg,1000)))),'r','linewidth',2)
% hold off
grid on
axis([-15 15 -100 10])
title('Spectrum of 4-Path Analysis Filter Bank Embedded In 30-Path Synthesis Filter Bank')
xlabel('Frequency')
ylabel('Log Mag (dB)')


figure(5)
dd=0.66;

subplot(4,1,1)
plot(0:519,5*real(exp(j*2*pi*dd)*x2(1:520)),'r','linewidth',2)

grid on
axis([0 450 -1.0 1.0])
title('Impulse Response (Real Part) Cascade Polyphase Filter Banks, 24 M&A per Output Sample')
ylabel('Amplitude')

subplot(4,1,2)
plot(0:298,5*real(gg),'linewidth',2)
grid on
axis([0 450 -1.0 1.0])
title('Impulse Response (Real Part) 299-Complex Tap FIR Filter, 478 M&A Per Output Sample')
ylabel('Amplitude')
xlabel('Time Index')

subplot(4,1,3)
plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(x2,1000)))),'r','linewidth',2)
% hold on
% plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(gg,1000)))),'linewidth',2)
% hold off
grid on
axis([-15 15 -100 10])
title('Frequency Response, Cascade Polyphase Filter Bank')
ylabel('Log Mag (dB)')

subplot(4,1,4)
plot((-0.5:1/1000:0.5-1/1000)*30,fftshift(20*log10(abs(fft(gg,1000)))),'linewidth',2)
grid on
axis([-15 15 -100 10])
title('Frequency Response, FIR Filter With Same Specifications')
ylabel('Log Mag (dB)')
xlabel('Frequency')
