% EE_656_HW_2b_2017, 8-path 8-to-1 down sampling filter, 16-taps per path
% Windowed sinc algorithm filter design, 0.1 dB Passband, 80 dB stopband

N1=(48/2)*80/22% 87
N1=1.4*N1  % 122... use 119... just missed, use 127
f1=1;
f2=3.0;
%h0=remez(86,[0 f1 f2 24]/24,{'myfrf',[1 1 0 0]},[1 120]);
phi=2*pi*(-63:63)*2/48;
h0=sin(phi)./phi;
h0(64)=1;
h0=h0.*kaiser(127,8)';
figure(21)
subplot(3,1,1)
plot(-63:63,h0,'linewidth',2);
grid on
axis([-65 65 -0.3 1.2])
title(['127-Tap Nyquist Filter, Windowed Sinc, 0.001-dB BW = \pm',num2str(f1),' kHZ, -80 dB BW =\pm',num2str(f2),' kHz, f_S = 48 kHz'],'fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

subplot(3,1,2)
fh0=fftshift(20*log10(abs(fft(h0*4/48,2048))));
plot((-0.5:1/2048:0.5-1/2048)*48,fh0,'linewidth',2)
hold on
plot( [-f1 -f1 +f1 +f1],[-90 -0.1 -0.1 -90],'r--','linewidth',2)
plot( [f2 f2 24],[-20 -80 -80],'r--','linewidth',2)
plot(-[f2 f2 24],[-20 -80 -80],'r--','linewidth',2)
hold off
grid on
axis([-24 24 -100 10])
title(['Spectrum; Passband 0-to-',num2str(f1),' kHz, Stopband ',num2str(f2),'-to-24 kHz, 80 db Attenuation'],'fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)

subplot(3,2,5)
plot((-0.5:1/2048:0.5-1/2048)*48,fh0,'linewidth',2)
hold on
plot( [-f1 -f1 +f1 +f1],[-0.0015 -0.001 -0.001 -0.0015],'r--','linewidth',2)
plot( [-f1 -f1 +f1 +f1],[+0.0015 +0.001 +0.001 +0.0015],'r--','linewidth',2)
hold off
grid on
axis([-f1-1 f1+1 -0.0015 0.0015])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)

subplot(3,2,6)
plot((-0.5:1/2048:0.5-1/2048)*48,fh0,'linewidth',2)
hold on
plot( [f1-1 +f1 +f1],[-0.1 -0.1 -90],'r--','linewidth',2)
plot( [f2 f2 20],[-20 -80 -80],'r--','linewidth',2)
hold off
grid on
axis([f1-1  f2+2 -100 5])
title('Transition BW and Stopband Detail','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)
text(1.1,-30,'0.001 dB BW','fontsize',14)
text(1.2,-45,[num2str(f1,'%5.1f'),' kHz'],'fontsize',14)
text(3.1,-70,['-80 dB BW = ',num2str(f2,'%5.1f'),' kHz'],'fontsize',14)


%%%%%%%%%%%%%%% 8-to-1 downsampling filter, 
% input sample rate = 48, output sample rate = 6
%%%%%%%%%%%%%%%%%%%%%%%% polyphase filter %%%%%%%%%%%%%%% 

hh0=reshape([h0 0],8,16);

reg=zeros(8,16);
v0=zeros(1,8)';
v1=zeros(8,34);
x0=zeros(1,200);
x0(1)=1;
x1=zeros(1,34);

m=1;
for n=1:8:200-8
    v0=fliplr(x0(n:n+7)).';
    reg=[v0 reg(:,1:15)];
        for k=1:8
          v1(k,m)=reg(k,:)*hh0(k,:)';
        end
       x1(m)=sum(v1(:,m));
        m=m+1;
end

figure(22)
subplot(3,1,1)
plot(0:15,x1(1:16),'-o','linewidth',2)
grid on
axis([0 13 -0.3 1.2])
title(['8-Path Filter, 8-to-1 Downsampled Impulse Response, 0.001-dB BW = \pm',num2str(f1),' kHZ, Transition BW = 2 kHz, Output sample rate f_S = 6 kHz'],'fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

subplot(3,1,2)
fx1=fftshift(20*log10(abs(fft(x1(1:33)/1.5,512))));  % sum(x1)=1.5001
plot((-0.5:1/512:0.5-1/512)*6,fx1,'linewidth',2)
hold on
plot([-1 -1 +1 +1],[-100 -0.001 -0.001 -100],'r--','linewidth',2)
hold off
grid on
axis([-3 3 -100 10])
title(['Spectrum; Passband \pm',num2str(f1),' kHz, Transition BW = 2 kHz'],'fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)

subplot(3,1,3)
plot((-0.5:1/512:0.5-1/512)*6,fx1-0.000,'linewidth',2)
hold on
plot([-1 -1 +1 +1],[-0.002 -0.001 -0.001 -0.002],'r--','linewidth',2)
plot([-1 -1 +1 +1],-[-0.002 -0.001 -0.001 -0.002],'r--','linewidth',2)
hold off
grid on
axis([-1.5 1.5 -0.002 0.002])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reg=zeros(8,16);
v0=zeros(1,8)';
v1=zeros(8,34);
x0=exp(j*2*pi*(0:1600)*1/256);
x1=zeros(1,200);

m=1;
for n=1:8:1600-8
    v0=fliplr(x0(n:n+7)).';
    reg=[v0 reg(:,1:15)];
        for k=1:8
          v1(k,m)=reg(k,:)*hh0(k,:)';
        end
       x1(m)=sum(v1(:,m));
        m=m+1;
end

figure(23)
scl=1/1.5;
for k=1:8
subplot(3,4,k)
plot(scl*real(v1(k,:)),'linewidth',2)
hold on
plot(scl*imag(v1(k,:)),'r','linewidth',2)
hold off
grid on
axis([0 100 -1.2 1.2])
title(['In-Band Tone Response Path(',num2str(k-1),')'],'fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)
end

subplot(3,3,7)
plot(scl*real(x1)/8,'linewidth',2)
hold on
plot(scl*imag(x1)/8,'r','linewidth',2)
hold off
grid on
axis([0 100 -1.2 1.2])
title('In-Band Tone Response Sum of 8-Paths','fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)

subplot(3,3,8)
ww=kaiser(128,10)';
ww=ww/sum(ww);
fx1=fftshift(20*log10(abs(fft(scl*x1(40:167)/1024,128))));
plot((-0.5:1/128:0.5-1/128)*6,fx1,'linewidth',2)
grid on
axis([-1 1 -100 10])
title('Frequency Response, Sum of 8-Paths','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)

subplot(3,3,9)
plot(0,0)
hold on
for k=1:8
    plot(fftshift((fft(scl*v1(k,(40:167))/128))),'-o','linewidth',2);
end
hold off
grid on
axis('square')
axis([-1.5 1.5 -1.5 1.5])
title('Nyquist Plot, Each of 8-Paths','fontsize',14)
ylabel('Real','fontsize',14)
xlabel('Imaginary','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reg=zeros(8,16);
v0=zeros(1,8)';
v1=zeros(8,34);
x0=exp(j*2*pi*(0:1600)*(6/48+1/256));
x1=zeros(1,200);

m=1;
for n=1:8:1600-8
    v0=fliplr(x0(n:n+7)).';
    reg=[v0 reg(:,1:15)];
        for k=1:8
          v1(k,m)=reg(k,:)*hh0(k,:)';
        end
       x1(m)=sum(v1(:,m));
        m=m+1;
end

figure(24)
for k=1:8
subplot(3,4,k)
plot(scl*real(v1(k,:)),'linewidth',2)
hold on
plot(scl*imag(v1(k,:)),'r','linewidth',2)
hold off
grid on
axis([0 100 -1.2 1.2])
title(['Stop-Band Tone Response Path(',num2str(k-1),')'],'fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)
end

subplot(3,3,7)
plot(scl*real(x1)/8,'linewidth',2)
hold on
plot(scl*imag(x1)/8,'r','linewidth',2)
hold off
grid on
axis([0 100 -1.2/10 1.2/10])
title('Stop-Band Tone Response Sum of 8-Paths','fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)
text(25,0.09,'Destructive Cancellation of','fontsize',14)
text(30,0.06,'Aliased Bandstop Tones','fontsize',14)

subplot(3,3,8)
ww=kaiser(128,10)';
ww=ww/sum(ww);
fx1=fftshift(20*log10(abs(fft(scl*x1(40:167)/200,128))));
plot((-0.5:1/128:0.5-1/128)*6,fx1,'linewidth',2)
grid on
axis([-3 3 -100 10])
title('Frequency Response, Sum of 8-Paths','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
text(-1.6,-14,'Destructive Cancellation of','fontsize',14)
text(-1.5,-26,'Aliased Bandstop Tones','fontsize',14)

subplot(3,3,9)
plot(0,0)
hold on
for k=1:8
    plot(fftshift(scl*(fft(v1(k,40:167)/128))),'-o','linewidth',2);
end
hold off
grid on
axis('square')
axis([-1.5 1.5 -1.5 1.5])
title('Nyquist Plot, Each of 8-Paths','fontsize',14)
ylabel('Real','fontsize',14)
xlabel('Imaginary','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(25)
subplot(3,3,1)
plot(0:126,h0,'.-','linewidth',2)
grid on
axis([0 128 -0.3 1.2])
title('Impulse Response, Prototype Filter','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

for k=1:8
    subplot(3,3,k+1)
    hh=zeros(1,127);
    hh(k:8:127)=h0(k:8:127);
    stem(0:126,hh,'marker','.','linewidth',2)
    grid on
    axis([0 128 -0.3 1.2])
title(['Impulse Response, Path(',num2str(k-1),')'],'fontsize',14)
    xlabel('Time Index','fontsize',14)
    ylabel('Amplitude','fontsize',14)
end

figure(26)
subplot(3,3,1)
plot((-0.5:1/1000:0.5-1/1000)*48,fftshift(20*log10(abs(fft(h0,1000)*4/48))),'linewidth',2)
grid on
axis([-24 24 -90 5])
title('Frequency Response, Prototype Filter','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)

for k=1:8
    subplot(3,3,k+1)
    hh=zeros(1,127);
    hh(k:8:127)=scl*h0(k:8:127);
    plot((-0.5:1/1000:0.5-1/1000)*48,fftshift(20*log10(abs(fft(hh,1000)))),'linewidth',2)
grid on
axis([-24 24 -90 5])
title(['Frequency Response, Path(',num2str(k-1),')'],'fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)
end


figure(27)
subplot(3,3,1)
fhh0=unwrap(angle(fftshift(fft(h0,1000))))/(2*pi);
plot((-0.5:1/1000:0.5-1/1000)*48,fhh0-fhh0(501),'r','linewidth',2)
grid on
axis([-24 24 -4 4])
title('Phase Response, Prototype Filter','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('\phi(\theta)/2\pi','fontsize',14)

for k=1:8
    subplot(3,3,k+1)
    hh=zeros(1,127);
    hh(k:8:127)=h0(k:8:127);
    fhh=unwrap(angle(fftshift(fft(hh,1000))))/(2*pi);
    plot((-0.5:1/1000:0.5-1/1000)*48,fhh-fhh(501),'linewidth',2)
    hold on
    plot((-0.5:1/1000:0.5-1/1000)*48,fhh0-fhh0(501),'r','linewidth',2)
    hold off
    grid on
    axis([-24 24 -32 32])
    title(['Phase Response, Path(',num2str(k-1),')'],'fontsize',14)
    xlabel('Frequency','fontsize',14)
    ylabel('\phi(\theta)/2\pi','fontsize',14)
end

figure(28)
plot(0,0)
hold on
for k=1:8
    hh=zeros(1,127);
    hh(k:8:127)=h0(k:8:127);
    fhh=unwrap(angle(fftshift(fft(hh,952))))/(2*pi);
    plot((-0.5:1/952:0.5-1/952)*48,fhh-fhh(477),'linewidth',2)
   
end
hold off
grid on
 title('Phase Response, all 8 Path Filters','fontsize',14)
 xlabel('Frequency','fontsize',14)
 ylabel('\phi(\theta)/2\pi','fontsize',14)
 axis([-24 24 -32 32])

figure(29)
plot(0,0)
hold on
vv=(0.5:-1/1000:-0.5+1/1000)*63;
for k=1:8
    hh=zeros(1,127);
    hh(k:8:127)=h0(k:8:127);
   %hh=fftshift(hh);
    fhh=unwrap(angle(fftshift(fft(hh,1000))))/(2*pi);
    plot((-0.5:1/1000:0.5-1/1000)*48,fhh-fhh(501)-vv,'linewidth',2)
end
hold off
grid on
title('Phase Response, all 8-Filters, Detrended Phase Profiles by Subtracting Constant Phase Slope tp Show Phase Differences','fontsize',14)
xlabel('Frequency','fontsize',14)
 ylabel('\phi(\theta)/2\pi','fontsize',14)
 axis([-24 24 -4 4])

figure(30)
subplot(1,2,1)
plot(0,0)
hold on
for k=1:8
fhh0=unwrap(angle(fftshift(fft(hh0(k,:),1000))))/(2*pi);
plot((-0.5:1/1000:0.5-1/1000)*6,fhh0-fhh0(501),'linewidth',2)
end
hold off
grid on
axis([0 3 -4 0])
title('Phase Profiles, all 8-Filters, Down-Sampled 8-to-1','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('\phi(\theta)/2\pi','fontsize',14)
text(1.25,-0.75,'Phase Slope is Time Delay','fontsize',14)

subplot(1,2,2)
plot(0,0)
hold on
for k=1:8
      
   [gd,w]=grpdelay(hh0(k,:),1,1000,'whole');
   gd=fftshift(gd);
    plot((-0.5:1/1000:0.5-1/1000)*6,gd,'linewidth',2)
end
plot([-1 -1 ],[6.75 8.25],'r--','linewidth',2)
plot([ 1  1], [6.75 8.25],'r--','linewidth',2)
hold off
title('Group Delay (Time Delay) Profiles, all 8-Filters, Down-Sampled 8-to-1','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Samples','fontsize',14)
grid on
axis([-3 3 6.75 8.25])
text(-0.5,8.1,'Passband','fontsize',14)
