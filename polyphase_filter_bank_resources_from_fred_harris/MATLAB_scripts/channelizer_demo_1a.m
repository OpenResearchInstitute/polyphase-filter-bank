% channelizer demo 1

% 30 channel analysis channelizer, 15-to-1 down sample

h=sinc(-(4-1/30):1/30:(4-1/30)).*kaiser(239,8)';
hx=h.*exp(j*2*pi*(-119:119)*1/30);
figure(1)
subplot(2,1,1)
plot(-(4-1/30):1/30:(4-1/30),h,'linewidth',2)
grid on
axis([-4 4 -0.3 1.2])
title('Impulse Response, 239-Tap Prototype Nyquist Filter for 30-Path Analysis Channelizer, 8-Taps per Path')
xlabel('Time Index')
ylabel('Amplitude')

subplot(2,1,2)
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(h/sum(h),2000)))),'linewidth',2)
hold on
plot((-0.5:1/2000:0.5-1/2000)*30,fftshift(20*log10(abs(fft(hx/sum(h),2000)))),'linewidth',2)
plot([+0.4 +0.6],[-1 -1]*6.02,'r','linewidth',2)
plot([0.5 0.5],[-12 -0],'r','linewidth',2)
hold off
grid on
axis([-3 3 -100 10])
title('Frequency Response, Prototype Nyquist Filter for 30-Path Analysis Channelizer')
xlabel('Frequency')
ylabel('Log Magnitude (dB)')



frq=[-5.91 -2.2 1 3.1 5.3];
x0=zeros(1,20000);
for k=1:5
  x0=x0+exp(j*2*pi*(0:20000-1)*frq(k)/30+j*2*pi*rand(1));
end

h2=reshape([h 0],30,8);
reg1=zeros(30,16);    % notice difference, changed from zeros(30,8) to zeros(30,16)

m=1;  % input index offset
s=1;  % output index

v1=zeros(1,30)';
v2=zeros(1,30)';
v3=zeros(1,30)';
v4=zeros(30,600);
flg=0;
for n=1:15:length(x0)-14
    v1(1:15)=flipud(x0(n:n+14).');
    v1(16:30)=v1(1:15);
    reg1=[v1 reg1(:,1:15)];
    
    for k=1:15
        v2(k)=reg1(k,1:2:16)*h2(k,:)';
        v2(k+15)=reg1(k,2:2:16)*h2(k+15,:)';
    end
    if flg==0
        flg=1;
    else
        flg=0;
        v2=[v2(16:30);v2(1:15)];
    end
    
    v3=ifft(v2);
    
    v4(:,s)=fftshift(v3);
    s=s+1;
    m=m+15;
end

figure(2)
for k=1:30
    subplot(5,6,k)
    plot(0:40,real(v4(k,1:41)))
    hold on
    plot(0:40,imag(v4(k,1:41)),'r')
    hold off
    grid on
    axis([2 40 -1.1 1.1])
     if rem(k,6)==1
        ylabel('Amplitude')
    end
    if k>23
        xlabel('Time Index')
    end
    text(0.25,16,['bin (',num2str(k-16),')'])
    text(30,1.3,['bin',num2str(k-16)])
end

figure(3)
ww=kaiser(500,8)';
ww=ww/sum(ww);
w2=kaiser(5000,10)';
w2=w2/sum(w2);
subplot(6,1,1)
plot((-0.5:1/5000:0.5-1/5000)*30,fftshift(20*log10(abs(fft(x0(1:5000).*w2)))))
hold on
for k=1:30
    plot([-0.82 -0.5 -0.18 0.18 0.5 0.82]+(-16+k),[-80 -6 0 0 -6 -80],'--r','linewidth',2)
end
 plot([-0.82 -0.5 -0.18 0.18 0.5 0.82]+15,[-80 -6 0 0 -6 -80],'--r','linewidth',2)
hold off
grid on
axis([-15 15 -90 10]);
title('Input Spectrum')
xlabel('Frequency')
ylabel('Log Mag (dB)')

for k=1:30
    subplot(6,6,k+6)
    plot(-0.5:1/500:0.5-1/500,fftshift(20*log10(abs(fft(v4(k,1:500).*ww)))))
    grid on
    axis([-0.5 0.5 -90 10])
    if rem(k,6)==1
        ylabel('Log Mag (dB)')
    end
    if k>23
        xlabel('Frequency')
    end
    text(0.25,16,['bin (',num2str(k-16),')'])
end



%%%%%%%%%%%%%%%% repeat channelizer but now with signals with bandwidth
% narrow band signals with BW = 1/2 sampled at 30 

g1=sinc(-5:1/5:5).*kaiser(51,8)'; %Shaping filter, 5-samples per symbol
g1=[0 0 g1 0 0];
g2=remez(70,[0 0.5 4.5 30]/30,{'myfrf',[1 1 0 0]},[1 10]); % 1-to-12 interpolating filter
g2=[0 6*g2];
x0=zeros(1,42000);
frq=[-6 -2 1 3 7];
for k=1:5
    xx0=(floor(2*rand(1,700))-0.5)/0.5+j*(floor(2*rand(1,700))-0.5)/0.5;
    reg_a=zeros(1,11);
    m=0;
    for n=1:700
        reg_a=[xx0(n) reg_a(1:10)];
        for nn=1:5
           xx1(m+nn)=reg_a*g1(nn:5:55)';
        end     
        m=m+5;
    end
    reg_b=zeros(1,6);
    m=0;
    for n=1:3500
        reg_b=[xx1(n) reg_b(1:5)];
        for nn=1:12
           xx2(m+nn)=reg_b*g2(nn:12:72)';
        end     
        m=m+12;
    end
x0=x0+xx2.*exp(j*2*pi*(1:42000)*frq(k)/30);
end


h2=reshape([h 0],30,8);
reg1=zeros(30,16);    % notice difference, changed from zeros(30,8) to zeros(30,16)

m=1;  % input index offset
s=1;  % output index

v1=zeros(1,30)';
v2=zeros(1,30)';
v3=zeros(1,30)';
v4=zeros(30,600);
flg=0;
for n=1:15:length(x0)-14
    v1(1:15)=flipud(x0(n:n+14).');
    v1(16:30)=v1(1:15);
    reg1=[v1 reg1(:,1:15)];
    
    for k=1:15
        v2(k)=reg1(k,1:2:16)*h2(k,:)';
        v2(k+15)=reg1(k,2:2:16)*h2(k+15,:)';
    end
    if flg==0
        flg=1;
    else
        flg=0;
        v2=[v2(16:30);v2(1:15)];
    end
    
    v3=ifft(v2);
    
    v4(:,s)=fftshift(v3);
    s=s+1;
    m=m+15;
end
% 
figure(4)
for k=1:30
    subplot(5,6,k)
    plot(0:200,real(v4(k,1:201)))
    hold on
    plot(0:200,imag(v4(k,1:201)),'r')
    hold off
    grid on
    axis([0 200 -1.1 1.1])
     if rem(k,6)==1
        ylabel('Amplitude')
    end
    if k>23
        xlabel('Time Index')
    end
    text(0.25,16,['bin (',num2str(k-16),')'])
    text(30,1.3,['bin',num2str(k-16)])
end

figure(5)
ww=kaiser(500,8)';
ww=10*ww/sum(ww);
w2=kaiser(5000,10)';
w2=10*w2/sum(w2);
subplot(6,1,1)
plot((-0.5:1/5000:0.5-1/5000)*30,fftshift(20*log10(abs(fft(x0(1:5000).*w2)))))
hold on
for k=1:30
    plot([-0.82 -0.5 -0.18 0.18 0.5 0.82]+(-16+k),[-80 -6 0 0 -6 -80],'--r','linewidth',2)
end
 plot([-0.82 -0.5 -0.18 0.18 0.5 0.82]+15,[-80 -6 0 0 -6 -80],'--r','linewidth',2)
hold off
grid on
axis([-15 15 -90 10]);
title('Input Spectrum')
xlabel('Frequency')
ylabel('Log Mag (dB)')

for k=1:30
    subplot(6,6,k+6)
    plot(-0.5:1/500:0.5-1/500,fftshift(20*log10(abs(fft(v4(k,1:500).*ww)))))
    grid on
    axis([-0.5 0.5 -90 10])
    if rem(k,6)==1
        ylabel('Log Mag (dB)')
    end
    if k>23
        xlabel('Frequency')
    end
    text(0.25,16,['bin (',num2str(k-16),')'])
end
% 

figure(7)
for k=1:30
    subplot(5,6,k)
    plot(v4(k,4:4:600),'r.')
        grid on
        axis('equal')
    axis([-1.1 1.1 -1.1 1.1])
     if rem(k,6)==1
        ylabel('Q Component')
    end
    if k>23
        xlabel('I Component')
    end
    text(0.1,1.3,['bin (',num2str(k-16),')'])
end

