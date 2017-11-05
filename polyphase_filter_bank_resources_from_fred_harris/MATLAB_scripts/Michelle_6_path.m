% N1=(36/2)*80/22 =65
h0=remez(64,[0 1.5 4.5 18]/18,{'myfrf',[1 1 0 0]},[1 1]);
f1=1.5;
f2=4.5;
% figure(1)
% subplot(3,1,1)
% plot(-32:32,h0/max(h0),'linewidth',2);
% grid on
% axis([-35 35 -0.3 1.2])
% title(['65-Tap Nyquist Filter, Windowed Sinc, 0.1-dB BW = \pm',num2str(f1),' kHZ, -80 dB BW =\pm',num2str(f2),' kHz, f_S = 36 kHz'],'fontsize',14)
% xlabel('Time Index','fontsize',14)
% ylabel('Amplitude','fontsize',14)
% text(15,0.7,['0.1 dB BW = ',num2str(f1,'%5.1f'),' kHz'],'fontsize',14)
% 
% subplot(3,1,2)
% fh0=fftshift(20*log10(abs(fft(h0,2048))));
% plot((-0.5:1/2048:0.5-1/2048)*36,fh0,'linewidth',2)
% hold on
% plot( [-f1 -f1 +f1 +f1],[-90 -0.1 -0.1 -90],'r--','linewidth',2)
% plot( [f2 f2 20],[-20 -80 -80],'r--','linewidth',2)
% plot(-[f2 f2 20],[-20 -80 -80],'r--','linewidth',2)
% hold off
% grid on
% axis([-18 18 -100 10])
% title(['Spectrum; Passband 0-to-',num2str(f1),' kHz, Stopband ',num2str(f2),'-to-18 kHz, 80 db Attenuation'],'fontsize',14)
% xlabel('Frequency (kHz)','fontsize',14)
% ylabel('Log Mag (dB)','fontsize',14)
% 
% subplot(3,2,5)
% plot((-0.5:1/2048:0.5-1/2048)*36,fh0,'linewidth',2)
% hold on
% plot( [-f1 -f1 +f1 +f1],[-0.15 -0.1 -0.1 -0.15],'r--','linewidth',2)
% plot( [-f1 -f1 +f1 +f1],[+0.15 +0.1 +0.1 +0.15],'r--','linewidth',2)
% hold off
% grid on
% axis([-f1-1 f1+1 -0.15 0.15])
% title('Zoom to Passband Ripple','fontsize',14)
% xlabel('Frequency (kHz)','fontsize',14)
% ylabel('Log Mag (dB)','fontsize',14)
% 
% subplot(3,2,6)
% plot((-0.5:1/2048:0.5-1/2048)*36,fh0,'linewidth',2)
% hold on
% plot( [f1-1 +f1 +f1],[-0.1 -0.1 -90],'r--','linewidth',2)
% plot( [f2 f2 20],[-20 -80 -80],'r--','linewidth',2)
% hold off
% grid on
% axis([f1-1  f2+2 -100 5])
% title('Transition BW and Stopband Detail','fontsize',14)
% xlabel('Frequency (kHz)','fontsize',14)
% ylabel('Log Mag (dB)','fontsize',14)
% 

%%%%%%%%%%%%%%% 6-to-1 downsampling filter, output sample rate = 6
%%%%%%%%%%%%%%%%%%%%%%%% polyphase filter %%%%%%%%%%%%%%% 
scl=max(h0);
hh0=reshape([h0/scl 0],6,11);
% 
% reg=zeros(6,11);
% v0=zeros(1,6)';
% v1=zeros(6,34);
% x0=zeros(1,200);
% x0(4)=1;
% x1=zeros(1,34);
% 
% m=1;
% for n=1:6:200-5
%     v0=fliplr(x0(n:n+5)).';
%     reg=[v0 reg(:,1:10)];
%         for k=1:6
%           v1(k,m)=reg(k,:)*hh0(k,:)';
%         end
%        x1(m)=sum(v1(:,m));
%         m=m+1;
% end
% 
% figure(2)
% subplot(3,1,1)
% plot(0:33,x1,'linewidth',2)
% grid on
% axis([0 10 -0.3 1.2])
% title(['Impulse Response, Output of 6-Path Filter, Downsampled 6-to-1: 0.1-dB BW = \pm',num2str(f1),' kHZ, -80 dB BW =\pm',num2str(f2),' kHz, f_S = 6 kHz'],'fontsize',14)
% xlabel('Time Index','fontsize',14)
% ylabel('Amplitude','fontsize',14)
% 
% subplot(3,1,2)
% fx1=fftshift(20*log10(abs(fft(x1(1:33)/sum(x1),512))));
% plot((-0.5:1/512:0.5-1/512)*6,fx1,'linewidth',2)
% hold on
% plot([-1 -1 +1 +1],[-100 -0.1 -0.1 -100],'r--','linewidth',2)
% hold off
% grid on
% axis([-3 3 -100 10])
% title(['Spectrum; Passband 0-to-',num2str(f1),' kHz, Transition Bandwidth ',num2str(f1),'-to-3 kHz, 80 db Attenuation'],'fontsize',14)
% xlabel('Frequency (kHz)','fontsize',14)
% ylabel('Log Mag (dB)','fontsize',14)
% 
% subplot(3,1,3)
% plot((-0.5:1/512:0.5-1/512)*6,fx1-0.103,'linewidth',2)
% hold on
% plot([-1 -1 +1 +1],[-0.2 -0.1 -0.1 -0.2],'r--','linewidth',2)
% plot([-1 -1 +1 +1],-[-0.2 -0.1 -0.1 -0.2],'r--','linewidth',2)
% hold off
% grid on
% axis([-3 3 -0.2 0.2])
% title('Zoom to Passband Ripple','fontsize',14)
% xlabel('Frequency (kHz)','fontsize',14)
% ylabel('Log Mag (dB)','fontsize',14)
% 

scl=max(h0);
hh0=reshape([h0/scl 0],6,11);
f_h0=fftshift(abs(fft(h0,128)));
w=kaiser(256,5')';

for k=0:1:256
reg=zeros(3,22);
v0=zeros(1,3)';
v1=zeros(1,6)';
v2=zeros(6,34);
%x0=exp(j*2*pi*(0:1200)*1/256)*exp(-j*2*pi*0.365);
x0=exp(j*2*pi*(0:1200)*k/256);

x1=zeros(1,200);

m=1;
flg=0;
for n=1:3:1200-3
    v0(1:3)=fliplr(x0(n:n+2)).';
    reg=[v0 reg(:,1:21)];
        for k=1:3
          v1(k)=reg(k,1:2:22)*hh0(k,:)';
          v1(k+3)=reg(k,2:2:22)*hh0(k+3,:)';
        end
        if flg==0
            flg=1;
        else
            flg=0;
            v1=[v1(4:6);v1(1:3)];
        end
        v2(:,m)=v1;
       x1(m)=sum(v2(:,m));
        m=m+1;
end

figure(3)
subplot(4,1,1)
plot((-0.5:1/512:0.5-1/512)*36,fftshift(abs(fft(x0(1:512)/512))),'linewidth',2.5)
hold on
for zz=-3:3
plot((-0.5:1/128:0.5-1/128)*36+zz*6,f_h0,'--','linewidth',2)
end
hold off
grid on
axis([-18 18 0 1.1])

for k=1:6
subplot(4,3,k+3)
plot(real(v2(k,:)),'linewidth',2.5)
hold on
plot(imag(v2(k,:)),'r','linewidth',2.5)
hold off
axis([0 100 -1.2 1.2])
grid on
title(['In-Band Tone Response Path(',num2str(k-1),')'],'fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)
end

subplot(4,3,10)
plot(scl*real(x1),'linewidth',2.5)
hold on
plot(scl*imag(x1),'r','linewidth',2.5)
hold off
grid on
axis([0 100 -1.2 1.2])
title('In-Band Tone Response Sum of 6-Paths','fontsize',14)
ylabel('Amplitude','fontsize',14)
xlabel('Time Index','fontsize',14)

subplot(4,3,11)
ww=kaiser(128,0)';
ww=ww/sum(ww);
fx1=fftshift(20*log10(abs(fft(scl*x1(21:148+128)/256,256))));
plot((-0.5:1/256:0.5-1/256)*6,fx1,'linewidth',2.5)
grid on
axis([-3 3 -120 10])
title('Frequency Response, Sum of 6-Paths','fontsize',14)
ylabel('Log Mag (dB)','fontsize',14)
xlabel('Frequency (kHz)','fontsize',14)

subplot(4,3,12)
plot(0,0)
hold on
for k=1:6
    plot(fftshift((fft(v2(k,21:148+128)/256))),'-o','linewidth',2.5);
end
hold off
grid on
axis('square')
axis([-1.5 1.5 -1.5 1.5])
title('Nyquist Plot, Each of 6-Paths','fontsize',14)
ylabel('Real','fontsize',14)
xlabel('Imaginary','fontsize',14)

pause(0.2)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% scl=max(h0);
% hh0=reshape([h0/scl 0],6,11);
% 
% reg=zeros(6,11);
% v0=zeros(1,6)';
% vv1=zeros(6,34);
% xx0=exp(j*2*pi*(0:1200)*(6/36+1/256))*exp(-j*2*pi*0.365);
% xx1=zeros(1,200);
% 
% m=1;
% for n=1:6:1200-5
%     v0=fliplr(xx0(n:n+5)).';
%     reg=[v0 reg(:,1:10)];
%         for k=1:6
%           vv1(k,m)=reg(k,:)*hh0(k,:)';
%         end
%        xx1(m)=sum(vv1(:,m));
%         m=m+1;
% end
% 
% figure(4)
% for k=1:6
% subplot(3,3,k)
% plot(real(vv1(k,:)),'linewidth',2.5)
% hold on
% plot(imag(vv1(k,:)),'r','linewidth',2.5)
% hold off
% grid on
% title(['Stop-Band Tone Response Path(',num2str(k-1),')'],'fontsize',14)
% ylabel('Amplitude','fontsize',14)
% xlabel('Time Index','fontsize',14)
% 
% end
% 
% subplot(3,3,7)
% plot(scl*real(xx1),'linewidth',2.5)
% hold on
% plot(scl*imag(xx1),'r','linewidth',2.5)
% hold off
% grid on
% axis([0 200 -1.2 1.2])
% title('Stop-Band Tone Response Sum of 6-Paths','fontsize',14)
% ylabel('Amplitude','fontsize',14)
% xlabel('Time Index','fontsize',14)
% text(45,0.9,'Destructive Cancellation of','fontsize',14)
% text(50,0.6,'Aliased Bandstop Tones','fontsize',14)
% 
% subplot(3,3,8)
% ww=kaiser(128,10)';
% ww=ww/sum(ww);
% fxx1=fftshift(20*log10(abs(fft(scl*xx1(21:148)/200,128))));
% plot((-0.5:1/128:0.5-1/128)*6,fxx1,'linewidth',2.5)
% grid on
% axis([-3 3 -100 10])
% title('Frequency Response, Sum of 6-Paths','fontsize',14)
% ylabel('Log Mag (dB)','fontsize',14)
% xlabel('Frequency (kHz)','fontsize',14)
% text(-1.6,-14,'Destructive Cancellation of','fontsize',14)
% text(-1.5,-26,'Aliased Bandstop Tones','fontsize',14)
% 
% subplot(3,3,9)
% plot(0,0)
% hold on
% for k=1:6
%     plot(fftshift((fft(vv1(k,21:148)/256))),'-o','linewidth',2.5);
% end
% hold off
% grid on
% axis('square')
% axis([-1.5 1.5 -1.5 1.5])
% title('Nyquist Plot, Each of 6-Paths','fontsize',14)
% ylabel('Real','fontsize',14)
% xlabel('Imaginary','fontsize',14)
% % 
% figure(5)
% subplot(3,2,1)
% plot((-0.5:1/768:0.5-1/768)*36,fftshift(20*log10(abs(fft(x0(1:768)/768)))),'linewidth',2)
% grid on
% axis([-8 8 -100 10])
% title('Spectrum, Input Signal, In-Band Tone','fontsize',14)
% xlabel('Frequency, (kHz)','fontsize',14)
% ylabel('Log Mag, (dB)','fontsize',14)
% 
% subplot(3,2,3)
% plot((-0.5:1/128:0.5-1/128)*6,fftshift(20*log10(abs(fft(v1(1,11:138)/768)))),'linewidth',2)
% grid on
% axis([-3 3 -100 10])
% title('Spectrum, Aliased Input Signal, In-Band Tone','fontsize',14)
% xlabel('Frequency, (kHz)','fontsize',14)
% ylabel('Log Mag, (dB)','fontsize',14)
% 
% subplot(3,2,5)
% plot((-0.5:1/128:0.5-1/128)*6,fftshift(20*log10(abs(fft(x1(11:138)/1024)))),'linewidth',2)
% grid on
% title('Spectrum, Output Signal, In-Band Tone','fontsize',14)
% xlabel('Frequency, (kHz)','fontsize',14)
% ylabel('Log Mag, (dB)','fontsize',14)
% axis([-3 3 -100 10])
% 
% subplot(3,2,2)
% plot((-0.5:1/768:0.5-1/768)*36,fftshift(20*log10(abs(fft(xx0(1:768)/768)))),'linewidth',2)
% grid on
% title('Spectrum, Input Signal, Out-of-Band Tone','fontsize',14)
% xlabel('Frequency, (kHz)','fontsize',14)
% ylabel('Log Mag, (dB)','fontsize',14)
% axis([-8 8 -100 10])
% 
% subplot(3,2,4)
% plot((-0.5:1/128:0.5-1/128)*6,fftshift(20*log10(abs(fft(vv1(1,11:138)/768)))),'linewidth',2)
% grid on
% axis([-3 3 -100 10])
% title('Spectrum, Aliased Input Signal, Out-of-Band Tone','fontsize',14)
% xlabel('Frequency, (kHz)','fontsize',14)
% ylabel('Log Mag, (dB)','fontsize',14)
% 
% subplot(3,2,6)
% plot((-0.5:1/128:0.5-1/128)*6,fftshift(20*log10(abs(fft(xx1(11:138)/1024)))),'linewidth',2)
% grid on
% title('Spectrum, Output Signal, Out-of-Band Tone','fontsize',14)
% xlabel('Frequency, (kHz)','fontsize',14)
% ylabel('Log Mag, (dB)','fontsize',14)
% axis([-3 3 -100 10])
% 
% 
% 
