% EE_656_HW_2b_2017, 8-path 8-to-1 down sampling filter, 11-taps per path
% REMEZ algorithm filter design, 0.1 dB Passbab, 80 dB stopband

% N1=(48/2)*80/22 = 87
f1=1;
f2=3.0;
h0=remez(86,[0 f1 f2 24]/24,{'myfrf',[1 1 0 0]},[1 120]);

figure(11)
subplot(3,1,1)
plot(-43:43,h0/max(h0),'linewidth',2);
grid on
axis([-50 50 -0.3 1.2])
title(['87-Tap Nyquist Filter, Windowed Sinc, 0.1-dB BW = \pm',num2str(f1),' kHZ, -80 dB BW =\pm',num2str(f2),' kHz, f_S = 48 kHz'],'fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

subplot(3,1,2)
fh0=fftshift(20*log10(abs(fft(h0,2048))));
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
plot( [-f1 -f1 +f1 +f1],[-0.15 -0.1 -0.1 -0.15],'r--','linewidth',2)
plot( [-f1 -f1 +f1 +f1],[+0.15 +0.1 +0.1 +0.15],'r--','linewidth',2)
hold off
grid on
axis([-f1-1 f1+1 -0.15 0.15])
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
text(1.1,-30,'0.1 dB BW','fontsize',14)
text(1.2,-45,[num2str(f1,'%5.1f'),' kHz'],'fontsize',14)
text(3.1,-70,['-80 dB BW = ',num2str(f2,'%5.1f'),' kHz'],'fontsize',14)


%%%%%%%%%%%%%%% 8-to-1 downsampling filter, 
% input sample rate = 48, output sample rate = 6
%%%%%%%%%%%%%%%%%%%%%%%% polyphase filter %%%%%%%%%%%%%%% 
scl=max(h0);
hh0=reshape([h0/scl 0],8,11);

reg=zeros(8,15);
v0=zeros(1,8)';
v1=zeros(8,34);
x0=zeros(1,200);
x0(5)=1;
x1=zeros(1,34);

m=1;
for n=1:8:200-8
    v0=fliplr(x0(n:n+7)).';
    reg=[v0 reg(:,1:10)];
        for k=1:8
          v1(k,m)=reg(k,:)*hh0(k,:)';
        end
       x1(m)=sum(v1(:,m));
        m=m+1;
end

figure(12)
subplot(3,1,1)
plot(0:10,x1(1:11),'-o','linewidth',2)
grid on
axis([0 10 -0.3 1.2])
title(['8-Path Filter, 8-to-1 Downsampled Impulse Response, 0.1-dB BW = \pm',num2str(f1),' kHZ, Transition BW = 2 kHz, Output sample rate f_S = 6 kHz'],'fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

subplot(3,1,2)
fx1=fftshift(20*log10(abs(fft(x1(1:33)/sum(x1),512))));
plot((-0.5:1/512:0.5-1/512)*6,fx1,'linewidth',2)
hold on
plot([-1 -1 +1 +1],[-100 -0.1 -0.1 -100],'r--','linewidth',2)
hold off
grid on
axis([-3 3 -100 10])
title(['Spectrum; Passband \pm',num2str(f1),' kHz, Transition BW = 2 kHz'],'fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)

subplot(3,1,3)
plot((-0.5:1/512:0.5-1/512)*6,fx1-0.1,'linewidth',2)
hold on
plot([-1 -1 +1 +1],[-0.2 -0.1 -0.1 -0.2],'r--','linewidth',2)
plot([-1 -1 +1 +1],-[-0.2 -0.1 -0.1 -0.2],'r--','linewidth',2)
hold off
grid on
axis([-3 3 -0.2 0.2])
title('Zoom to Passband Ripple','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)


% scl=max(h0);
% hh0=reshape([h0/scl 0],8,11);

reg=zeros(8,11);
v0=zeros(1,8)';
v1=zeros(8,34);
x0=exp(j*2*pi*(0:1600)*1/256)*exp(-j*2*pi*0.3125);
x1=zeros(1,200);

m=1;
for n=1:8:1600-8
    v0=fliplr(x0(n:n+7)).';
    reg=[v0 reg(:,1:10)];
        for k=1:8
          v1(k,m)=reg(k,:)*hh0(k,:)';
        end
       x1(m)=sum(v1(:,m));
        m=m+1;
end

figure(13)
for k=1:8
subplot(3,4,k)
plot(real(v1(k,:)),'linewidth',2)
hold on
plot(imag(v1(k,:)),'r','linewidth',2)
hold off
grid on
title(['In-Band Tone Response Path(',num2str(k-1),')'],'fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)
end

subplot(3,3,7)
plot(scl*real(x1),'linewidth',2)
hold on
plot(scl*imag(x1),'r','linewidth',2)
hold off
grid on
axis([0 200 -1.2 1.2])
title('In-Band Tone Response Sum of 8-Paths','fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)

subplot(3,3,8)
ww=kaiser(128,10)';
ww=ww/sum(ww);
fx1=fftshift(20*log10(abs(fft(scl*x1(21:148)/200,128))));
plot((-0.5:1/128:0.5-1/128)*6,fx1,'linewidth',2)
grid on
axis([-3 3 -100 10])
title('Frequency Response, Sum of 8-Paths','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)

subplot(3,3,9)
plot(0,0)
hold on
for k=1:8
    plot(fftshift((fft(v1(k,21:148)/256))),'-o','linewidth',2);
end
hold off
grid on
axis('square')
axis([-1.5 1.5 -1.5 1.5])
title('Nyquist Plot, Each of 8-Paths','fontsize',14)
ylabel('Real','fontsize',14)
xlabel('Imaginary','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% scl=max(h0);
% hh0=reshape([h0/scl 0],8,15);

reg=zeros(8,11);
v0=zeros(1,8)';
v1=zeros(8,34);
x0=exp(j*2*pi*(0:1600)*(6/48+1/256))*exp(-j*2*pi*0.3125);
x1=zeros(1,200);

m=1;
for n=1:8:1600-8
    v0=fliplr(x0(n:n+7)).';
    reg=[v0 reg(:,1:10)];
        for k=1:8
          v1(k,m)=reg(k,:)*hh0(k,:)';
        end
       x1(m)=sum(v1(:,m));
        m=m+1;
end

figure(14)
for k=1:8
subplot(3,4,k)
plot(real(v1(k,:)),'linewidth',2)
hold on
plot(imag(v1(k,:)),'r','linewidth',2)
hold off
grid on
title(['Stop-Band Tone Response Path(',num2str(k-1),')'],'fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)
end

subplot(3,3,7)
plot(scl*real(x1),'linewidth',2)
hold on
plot(scl*imag(x1),'r','linewidth',2)
hold off
grid on
axis([0 200 -1.2 1.2])
title('Stop-Band Tone Response Sum of 8-Paths','fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)
text(45,0.9,'Destructive Cancellation of','fontsize',14)
text(50,0.6,'Aliased Bandstop Tones','fontsize',14)

subplot(3,3,8)
ww=kaiser(128,10)';
ww=ww/sum(ww);
fx1=fftshift(20*log10(abs(fft(scl*x1(21:148)/200,128))));
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
    plot(fftshift((fft(v1(k,21:148)/256))),'-o','linewidth',2);
end
hold off
grid on
axis('square')
axis([-1.5 1.5 -1.5 1.5])
title('Nyquist Plot, Each of 8-Paths','fontsize',14)
ylabel('Real','fontsize',14)
xlabel('Imaginary','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(15)
subplot(3,3,1)
plot(0:86,h0,'.-','linewidth',2)
grid on
title('Impulse Response, Prototype Filter','fontsize',14)
xlabel('Time Index','fontsize',14)
ylabel('Amplitude','fontsize',14)

for k=1:8
    subplot(3,3,k+1)
    hh=zeros(1,87);
    hh(k:8:87)=h0(k:8:87);
    stem(0:86,hh,'marker','.','linewidth',2)
    grid on
    title(['Impulse Response, Path(',num2str(k-1),')'],'fontsize',14)
    xlabel('Time Index','fontsize',14)
    ylabel('Amplitude','fontsize',14)
end

figure(16)
subplot(3,3,1)
plot((-0.5:1/1000:0.5-1/1000)*48,fftshift(20*log10(abs(fft(h0,1000)))),'linewidth',2)
grid on
axis([-24 24 -90 5])
title('Frequency Response, Prototype Filter','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)

for k=1:8
    subplot(3,3,k+1)
    hh=zeros(1,87);
    hh(k:8:87)=8*h0(k:8:87);
    plot((-0.5:1/1000:0.5-1/1000)*48,fftshift(20*log10(abs(fft(hh,1000)))),'linewidth',2)
grid on
axis([-24 24 -90 5])
title(['Frequency Response, Path(',num2str(k-1),')'],'fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)

end


figure(17)
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
    hh=zeros(1,87);
    hh(k:8:87)=h0(k:8:87);
    fhh=unwrap(angle(fftshift(fft(hh,1000))))/(2*pi);
    plot((-0.5:1/1000:0.5-1/1000)*48,fhh-fhh(501),'linewidth',2)
    hold on
    plot((-0.5:1/1000:0.5-1/1000)*48,fhh0-fhh0(501),'r','linewidth',2)
    hold off
    grid on
    axis([-24 24 -24 24])
    title(['Phase Response, Path(',num2str(k-1),')'],'fontsize',14)
    xlabel('Frequency','fontsize',14)
    ylabel('\phi(\theta)/2\pi','fontsize',14)
end

figure(18)
plot(0,0)
hold on
for k=1:8
    hh=zeros(1,87);
    hh(k:8:87)=h0(k:8:87);
    fhh=unwrap(angle(fftshift(fft(hh,952))))/(2*pi);
    plot((-0.5:1/952:0.5-1/952)*48,fhh-fhh(477),'linewidth',2)
   
end
hold off
grid on
 title('Phase Response, all 8 Path Filters','fontsize',14)
 xlabel('Frequency','fontsize',14)
 ylabel('\phi(\theta)/2\pi','fontsize',14)
 axis([-24 24 -24 24])

figure(19)
plot(0,0)
hold on
vv=(0.5:-1/1000:-0.5+1/1000)*43;
for k=1:8
    hh=zeros(1,87);
    hh(k:8:87)=h0(k:8:87);
   %hh=fftshift(hh);
    fhh=unwrap(angle(fftshift(fft(hh,1000))))/(2*pi);
    plot((-0.5:1/1000:0.5-1/1000)*48,fhh-fhh(501)-vv,'linewidth',2)
end
hold off
grid on
title('Phase Response, all 8-Filters, Detrended Phase Profiles by Subtracting Constant Phase Slope tp Show Phase Differences','fontsize',14)
xlabel('Frequency','fontsize',14)
 ylabel('\phi(\theta)/2\pi','fontsize',14)
 axis([-24 24 -2.5 2.5])

figure(20)
subplot(1,2,1)
plot(0,0)
hold on
for k=1:8
fhh0=unwrap(angle(fftshift(fft(hh0(k,:),1000))))/(2*pi);
plot((-0.5:1/1000:0.5-1/1000)*6,fhh0-fhh0(501),'linewidth',2)
end
hold off
grid on
axis([-3 3 -3 3])
title('Phase Profiles, all 8-Filters, Down-Sampled 8-to-1','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('\phi(\theta)/2\pi','fontsize',14)

subplot(1,2,2)
plot(0,0)
hold on
for k=1:8
      
   [gd,w]=grpdelay(hh0(k,:),1,1000,'whole');
   gd=fftshift(gd);
    plot((-0.5:1/1000:0.5-1/1000)*6,gd,'linewidth',2)
end
plot([-1 -1 ],[4.4 5.5],'r--','linewidth',2)
plot([ 1  1], [4.4 5.5],'r--','linewidth',2)
hold off
title('Group Delay (Time Delay) Profiles, all 8-Filters, Down-Sampled 8-to-1','fontsize',14)
xlabel('Frequency','fontsize',14)
ylabel('Samples','fontsize',14)
grid on
axis([-3 3 4 6])
