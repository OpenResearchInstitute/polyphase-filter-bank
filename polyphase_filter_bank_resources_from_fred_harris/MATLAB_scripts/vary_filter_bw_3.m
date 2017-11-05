
clear all
vv=remez(54,[0 3.25 4 4 4.75 12]/12,{'myfrf',[1 1 0.50 0.50 0 0]}, [1 1 10]);

 
figure(4)
 alpha=-0.1;
 aa1=[1 -alpha];
 bb1=[-alpha 1];
 
 hh=1;
 rr=[1 zeros(1,100)];
 yy=zeros(1,101);
 
 for nn=1:54
     yy=yy+rr*vv(nn);
     rr=filter(bb1,aa1,rr);
 end
 yy=yy+rr*vv(55);

 subplot(3,3,1)
 plot(yy)
 grid on
 axis([0 60 -0.15 0.4]) 
title('Impulse Response, \alpha=-0.05')

subplot(3,3,4)
  plot((-0.5:1/1024:0.5-1/1024)*24,fftshift(20*log10(abs(fft(yy,1024)))))
 grid on
 axis([-12 12 -90 5])
title('Spectrum, Tuned to (\pm 4.0 MHz)')
xlabel('Frequency (MHz)')
ylabel('Log Magnitude (dB)')

 subplot(3,3,7)
  phi=unwrap(angle(fftshift(fft(yy,1024))));
  
  dphi=conv(phi,[-1 0 1])*512/(2*pi);
  dphi=[dphi(2:1024) dphi(1024)];
  
  m_yy=zeros(1,1024);
  m_yy((-200:200)+513)=1;
  plot((-0.5:1/1024:0.5-1/1024)*24,dphi.*m_yy)
  grid on
  axis([-12 12 20 35])
 title('Group Delay: 1 - Sample = 44.1 nsec')
  xlabel('Frequency (MHz)')
  ylabel('Delay (samples)')
  
 
  %%%%%%%%%%%%%%%%%%%%%
  alpha=-0.00;
 aa1=[1 -alpha];
 bb1=[-alpha 1];
 
 hh=1;
 rr=[1 zeros(1,100)];
 yy=zeros(1,101);
 
 for nn=1:54
     yy=yy+rr*vv(nn);
     rr=filter(bb1,aa1,rr);
 end
 yy=yy+rr*vv(55);

 subplot(3,3,2)
 plot(yy)
 grid on
 axis([0 60 -0.15 0.4]) 
title('Impulse Response, \alpha=0.00')

subplot(3,3,5)
  plot((-0.5:1/1024:0.5-1/1024)*24,fftshift(20*log10(abs(fft(yy,1024)))))
 grid on
 axis([-12 12 -90 5])
title('Spectrum, prototype (\pm 3.5 MHz)')
xlabel('Frequency (MHz)')
ylabel('Log Magnitude (dB)')

 subplot(3,3,8)
  phi=unwrap(angle(fftshift(fft(yy,1024))));
  
  dphi=conv(phi,[-1 0 1])*512/(2*pi);
  dphi=[dphi(2:1024) dphi(1024)];
  
  m_yy=zeros(1,1024);
  m_yy((-200:200)+513)=1;
  plot((-0.5:1/1024:0.5-1/1024)*24,dphi.*m_yy)
  grid on
  axis([-12 12 20 35])
 title('Group Delay: 1 - Sample = 44.1 nsec')
  xlabel('Frequency (MHz)')
  ylabel('Delay (samples)')
  
  
  %%%%%%%%%%%%%%%%%%%%%
  alpha=+0.1;
 aa1=[1 -alpha];
 bb1=[-alpha 1];
 
 hh=1;
 rr=[1 zeros(1,100)];
 yy=zeros(1,101);
 
 for nn=1:54
     yy=yy+rr*vv(nn);
     rr=filter(bb1,aa1,rr);
 end
 yy=yy+rr*vv(55);

 subplot(3,3,3)
 plot(yy)
 grid on
 axis([0 60 -0.15 0.4]) 
title('Impulse Response, \alpha=-0.05')

subplot(3,3,6)
  plot((-0.5:1/1024:0.5-1/1024)*24,fftshift(20*log10(abs(fft(yy,1024)))))
 grid on
 axis([-12 12 -90 5])
title('Spectrum, Tuned to (\pm 3.0 MHz)')
xlabel('Frequency (MHz)')
ylabel('Log Magnitude (dB)')

 subplot(3,3,9)
  phi=unwrap(angle(fftshift(fft(yy,1024))));
  
  dphi=conv(phi,[-1 0 1])*512/(2*pi);
  dphi=[dphi(2:1024) dphi(1024)];
  
  m_yy=zeros(1,1024);
  m_yy((-170:170)+513)=1;
  plot((-0.5:1/1024:0.5-1/1024)*24,dphi.*m_yy)
  grid on
  axis([-12 12 20 35])
 title('Group Delay: 1 - Sample = 44.1 nsec')
  xlabel('Frequency (MHz)')
  ylabel('Delay (samples)')
  
 