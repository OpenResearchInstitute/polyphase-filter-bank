% high quality interpolator
h0=remez(639,[0 0.25 0.75 32]/32,{'myfrf',[1 1 0 0]},[1,100]);
h0=h0*64;
h1x=conv(h0,[1 0 -1]/2);
h1=h1x(2:641);
h2x=conv(h1,[1 0 -1]/2);
h2=h2x(2:641);

g0=reshape(h0,64,10);
g1=reshape(h1,64,10);
g2=reshape(h2,64,10);

% low pass filter
hh=remez(150,[0 0.2 0.23 0.5]/0.5,{'myfrf',[1 1 0 0]},[1 400]);
fhh=fftshift(20*log10(abs(fft(hh,1024))));

x0=0.5*ones(1,2001);
for kk=1:8
    x0=x0+cos(2*pi*(0:2000)*kk/40+2*pi*rand(1));
end
%x0=cos(2*pi*(0:2000)*0.10);
%x0=[50*sinc(-50:.1:50).*kaiser(1001,10)' zeros(1,1000)];

vv=[0:.02:1 .98:-0.02:0.02 0:.02:1 .98:-0.02:0.02 0:.02:1];
    ww=kaiser(1024,15)';
    ww=ww/sum(ww);
for kk=1:251
incr1=64*(1+vv(kk)/1.8);
incr2=64/(1+vv(kk)/1.8);

reg1=zeros(1,10);
accum1=1;
mm=1;
nn=1;
flg1=1;
%y0_sv=zeros(1,4000);
%y1_sv=zeros(1,4000);
y2_sv=zeros(1,4000);
while nn<length(x0)-5;
    if flg1==1;
    reg1=[x0(nn) reg1(1:9)];
    nn=nn+1;
    flg1=0;
    end
    if accum1<65
    pntr1=floor(accum1);
    del1=accum1-pntr1;
    y0=reg1*g0(pntr1,:)';
    y1=reg1*g1(pntr1,:)';
    y2=reg1*g2(pntr1,:)';
    %y0_sv(mm)=y0;
    %y1_sv(mm)=y0+del1*y1;
    y2_sv(mm)=y0+del1*(y1+del1*y2/2);
    mm=mm+1;
    
    accum1=accum1+incr1;
    else
    if accum1>=65;
        accum1=accum1-64;
        flg1=1;
    end
    end
end
    y3=filter(hh,1,y2_sv);
reg1=zeros(1,10);
accum1=1;
mm=1;
nn=1;
flg1=1;
%y0_sv=zeros(1,4000);
%y1_sv=zeros(1,4000);
y4=zeros(1,4000);
while nn<length(y3)-5;
    if flg1==1;
    reg1=[y3(nn) reg1(1:9)];
    nn=nn+1;
    flg1=0;
    end
    if accum1<65
    pntr1=floor(accum1);
    del1=accum1-pntr1;
    y0=reg1*g0(pntr1,:)';
    y1=reg1*g1(pntr1,:)';
    y2=reg1*g2(pntr1,:)';
    %y0_sv(mm)=y0;
    %y1_sv(mm)=y0+del1*y1;
    y4(mm)=y0+del1*(y1+del1*y2/2);
    mm=mm+1;
    
    accum1=accum1+incr2;
    else
    if accum1>=65;
        accum1=accum1-64;
        flg1=1;
    end
    end
end    
    figure(2)
    subplot(4,1,1)
    plot(-0.5:1/1024:0.5-1/1024,fftshift(20*log10(abs(fft(x0(201:1224).*ww)))))
    hold on
    plot(-0.5:1/1024:0.5-1/1024,fhh,'r')
    hold off
    grid on
    axis([-0.5 0.5 -140 10])
    ylabel('Log Mag (dB)')
    title('Spectra, Input and Output')

    subplot(4,1,2)
    plot(-0.5:1/1024:0.5-1/1024,fftshift(20*log10(abs(fft(y2_sv(201:1224).*ww)))))
    hold on
    plot(-0.5:1/1024:0.5-1/1024,fhh,'r')
    hold off
    grid on
    axis([-0.5 0.5 -140 10])
    ylabel('Log Mag (dB)')
    
    subplot(4,1,3)
    plot(-0.5:1/1024:0.5-1/1024,fftshift(20*log10(abs(fft(y3(201:1224).*ww)))))
    grid on
    axis([-0.5 0.5 -140 10])
    ylabel('Log Mag (dB)')

    subplot(4,1,4)
    plot(-0.5:1/1024:0.5-1/1024,fftshift(20*log10(abs(fft(y4(201:1224).*ww)))))
    grid on
    axis([-0.5 0.5 -140 10])
    ylabel('Log Mag (dB)')
    pause(0.1)
end