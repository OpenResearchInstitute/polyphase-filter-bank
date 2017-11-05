
h=remez(98,[0 0.4 0.6 2.5]/2.5,[1 1 0 0],[1 100]);

hh=[h 0];
hh2=reshape(hh,5,20);

reg=zeros(5,20);

x0=zeros(1,12000);

ff=[0.05 0.15 0.25 0.35 0.40 0.42 0.44 0.46 0.48 0.51 0.53 0.55 0.57 0.59 0.61 0.65 0.7 1.01 1.49 1.8];
for k=1:19
    x0=x0+cos(2*pi*(1:12000)*(ff(k)/5) +2*pi*rand(1));
end

m=1;
for n=1:5:12000-4
    
    v1=fliplr(x0(n:n+4)).';
    reg=[v1 reg(:,1:19)];
    
    for k=1:5
        v2(k)=reg(k,:)*hh2(k,:)';
    end
    y(m)=sum(v2);
    m=m+1;
end

figure(1)
ww=kaiser(2000,10)';
ww=ww/sum(ww);
subplot(2,1,1)
plot((-0.5:1/2000:0.5-1/2000)*5,fftshift(20*log10(abs(fft(x0(1:2000).*ww)))),'linewidth',2)
hold on
plot((-0.5:1/2000:0.5-1/2000)*5,fftshift(20*log10(abs(fft(h,2000)))),'r','linewidth',2)

hold off
grid on
axis([-2.5 2.5 -100 10])

subplot(2,1,2)
plot((-0.5:1/2000:0.5-1/2000)*1,fftshift(20*log10(abs(fft(y(101:2100).*ww,2000)))),'linewidth',2)
grid on
axis([-0.5 0.5 -100 10])

