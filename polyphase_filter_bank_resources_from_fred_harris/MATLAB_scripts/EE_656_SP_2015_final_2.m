% EE_656_Spring_2015_Final_2

% multiple 5=MHz qpsk signals,5-MHz Symbol rates, 6-MHz Bandwidths, 80-MHz Sample Rate

%h1=nyq_2zz(1,16,0.2,18,1,0);
h1=sinc(-10+1/16:1/16:10-1/16).*kaiser(319,6)';
figure(1)
subplot(2,1,1)
plot(h1/max(h1))
grid on
axis([0 350 -0.3 1.2])

subplot(2,1,2)
plot((-0.5:1/4096:0.5-1/4096)*80,fftshift(20*log10(abs(fft(h1/sum(h1),4096)))))
hold on
plot([-8 -3 -3],[-60 -60 -20],'r')
plot([+8 +3 +3],[-60 -60 -20],'r')
hold off
grid on
axis([-8 8 -90 10])

x2=zeros(1,17600);

ff=[-3 -2 -1 0 1 2 3]/8;
aa=[ 1  1  1 1 1 0 1]/8;
%aa=[ 0 0 0 0 1 0 0 0]/8;

for k=1:7
    x0=(floor(2*rand(1,1100))-0.5)/0.5+j*(floor(2*rand(1,1100))-0.5)/0.5;
    x1=zeros(1,17600);
    x1(1:16:17600)=x0;
    x1=filter(h1,1,x1);
    x2=x2+aa(k)*(x1.*exp(j*2*pi*(0:17599)*ff(k)));
end

figure(2)
subplot(2,1,1)
plot(real(x2(1:1000)));
grid on
axis([0 500 -1.2 1.2])

ww=kaiser(4096,10)';
ww=100*ww/sum(ww);
w8=kaiser(2048,10)';
w8=100*w8/sum(w8);

subplot(2,1,2)
plot((-0.5:1/4096:0.5-1/4096)*80,fftshift(20*log10(abs(fft(x2(201:4296).*ww)))))
grid on
axis([-40 40 -90 10])

% Part 1. Separate Channel down conversion, filter and 8-to-1 downsample
g1=remez(60,[0 3 7 40]/40,{'myfrf',[1 1 0 0]},[1 10]);
x3=x2.*exp(-j*2*pi*(0:17599)*1/8);
x3a=filter(g1,1,x3);
x3b=x3a(14:8:end);

figure(3)
subplot(3,1,1)
plot(0:2000,8*real(x3a(1:2001)))
hold on
plot((0:16:2001)+12,8*real(x3b(1:2:251)),'ro')
hold off
grid on
axis([0 2000 -1.5 1.5])

subplot(3,1,2)
plot((-0.5:1/4096:0.5-1/4096)*80,fftshift(20*log10(abs(fft(x3a(201:4296).*ww)))))
hold on
plot((-0.5:1/2048:0.5-1/2048)*80,fftshift(20*log10(abs(fft(g1,2048)))),'r')
hold off
grid on
axis([-40 40 -90 10])

subplot(3,1,3)
plot((-0.5:1/2048:0.5-1/2048)*10,fftshift(20*log10(abs(fft(x3b(21:2068).*w8)))))
grid on
axis([-5 5 -90 10])


figure(4)
subplot(2,1,1)
plot(8*x3a(14:16:end),'r.')
axis('equal')
axis([-1.2 1.2 -1.2 1.2])
grid on
subplot(2,1,2)
plot(8*real(x3a(14:16:end)),'ro')
grid on
axis([0 500 -1.2 1.2])

% part 2, 10 path, polyphase filter, 10-to-1 Down sample 80-to-8 MHz
%         interpolator 8-to-10 MHz

%polyphase 10 path filter 
h2=remez(62,[0 3 7 40]/40,{'myfrf',[1 1 0 0]},[1 10]);
hh2=8*reshape([0 h2],8,8);
reg2=zeros(8,8);

v1=zeros(1,8)';
v2=zeros(1,8)';
v3=zeros(1,8)';
v4=zeros(8,2125);
reg=zeros(8,8);
m=1;
for n=1:8:17600-8
    v1=flipud(x2(n:n+7).');
    reg2=[v1 reg2(:,1:7)];
    for k=1:8
        v2(k)=reg2(k,:)*hh2(k,:)';
    end
    v3=8*ifft(fftshift(v2));
    v4(:,m)=v3;
    m=m+1;
end
% 
figure(5)
subplot(2,1,1)
plot(real(v4(2,1:500)))
grid on
axis([0 500 -2 2])

subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048)*10,fftshift(20*log10(abs(fft(v4(2,21:2068).*w8/4)))))
grid on
axis([-5 5 -90 10])


%%%%%%%%%%%%%% option 3 recursive Linear phase filter
%   24 Poles at Origin
%  
%   branch 1 - type1     coefficient    0.2869013219880045
%   branch 1 - type2 1st coefficient    -0.1785648602970922
%   branch 1 - type2 2nd coefficient    0.0226004172326217
% 			
%   branch 2 - type1     coefficient    0.4076140351608992
%   branch 2 - type2 1st coefficient    -0.1872725963096948
%   branch 2 - type2 2nd coefficient    0.0259854988225092
% 			
%   branch 3 - type1     coefficient    0.5111759052441421
%   branch 3 - type2 1st coefficient    -0.1747071895685166
%   branch 3 - type2 2nd coefficient    0.02462786753586059
% 			
%   branch 4 - type1     coefficient    0.6083890392634814
%   branch 4 - type2 1st coefficient    -0.1511534587793413
%   branch 4 - type2 2nd coefficient    0.02092667255651533
% 			
%   branch 5 - type1     coefficient    0.703646358874289
%   branch 5 - type2 1st coefficient    -0.1204118725193743
%   branch 5 - type2 2nd coefficient    0.01602156653771641
% 			
%   branch 6 - type1     coefficient    0.7994926311817401
%   branch 6 - type2 1st coefficient    -0.08434530463787691
%   branch 6 - type2 2nd coefficient    0.01060162856030615
% 			
%   branch 7 - type1     coefficient    0.897761894874325
%   branch 7 - type2 1st coefficient    -0.04399879293330994
%   branch 7 - type2 2nd coefficient    0.005141808905950611
  
coef3(1,:)=[0.2869013219880045   -0.1785648602970922   0.0226004172326217];
coef3(2,:)=[0.4076140351608992   -0.1872725963096948   0.0259854988225092];
coef3(3,:)=[0.5111759052441421   -0.1747071895685166   0.02462786753586059];
coef3(4,:)=[0.6083890392634814   -0.1511534587793413   0.02092667255651533];
coef3(5,:)=[0.703646358874289    -0.1204118725193743   0.01602156653771641];
coef3(6,:)=[0.7994926311817401   -0.08434530463787691  0.01060162856030615];
coef3(7,:)=[0.897761894874325    -0.04399879293330994  0.005141808905950611];

zz=zeros(1,7);

d0=[1 zeros(1,24)];
d11=[1 zz coef3(1,1)];
d12=[1 zz coef3(1,2) zz coef3(1,3)];
d21=[1 zz coef3(2,1)];
d22=[1 zz coef3(2,2) zz coef3(2,3)];
d31=[1 zz coef3(3,1)];
d32=[1 zz coef3(3,2) zz coef3(3,3)];
d41=[1 zz coef3(4,1)];
d42=[1 zz coef3(4,2) zz coef3(4,3)];
d51=[1 zz coef3(5,1)];
d52=[1 zz coef3(5,2) zz coef3(5,3)];
d61=[1 zz coef3(6,1)];
d62=[1 zz coef3(6,2) zz coef3(6,3)];
d71=[1 zz coef3(7,1)];
d72=[1 zz coef3(7,2) zz coef3(7,3)];

x0=[1 zeros(1,500)];
%x0=exp(j*2*pi*(1:17000)*1./8);


    y0=filter(fliplr(d0),d0,x0);
    y1=filter(fliplr([1 0]),[1 0],x0);
    y1=filter(fliplr(d11),d11,y1);
    y1=filter(fliplr(d12),d12,y1);
    
    y2=filter(fliplr([1 0 0]),[1 0 0],x0);
    y2=filter(fliplr(d21),d21,y2);
    y2=filter(fliplr(d22),d22,y2);
    
    y3=filter(fliplr([1 0 0 0]),[1 0 0 0],x0);
    y3=filter(fliplr(d31),d31,y3);
    y3=filter(fliplr(d32),d32,y3);
    
    y4=filter(fliplr([1 0 0 0 0]),[1 0 0 0 0],x0);
    y4=filter(fliplr(d41),d41,y4);
    y4=filter(fliplr(d42),d42,y4);
    
    y5=filter(fliplr([1 0 0 0 0 0]),[1 0 0 0 0 0],x0);
    y5=filter(fliplr(d51),d51,y5);
    y5=filter(fliplr(d52),d52,y5);
    
    y6=filter(fliplr([1 0 0 0 0 0 0]),[1 0 0 0 0 0 0],x0);
    y6=filter(fliplr(d61),d61,y6);
    y6=filter(fliplr(d62),d62,y6);
    
    y7=filter(fliplr([1 0 0 0 0 0 0 0]),[1 0 0 0 0 0 0 0],x0);
    y7=filter(fliplr(d71),d71,y7);
    y7=filter(fliplr(d72),d72,y7);
        
    h8=[y0+y1+y2+y3+y4+y5+y6+y7]/8;
    
x3=x2.*exp(-j*2*pi*(0:17599)*1/8);
x8a=filter(h8,1,x3);


% 10 path delay line
gg=remez(119,[0 4 8 64]/64,[1 1 0 0]);
dd=10*reshape(gg,10,12);

x8b=conv(x8a(14:8:end),dd(3,:));

 
figure(6)
subplot(3,1,1)
plot(0:200,h8(1:201))
grid on
axis([0 120 -0.05 0.15])

subplot(3,1,2)
fh8=fftshift(fft(h8,2048));
plot((-0.5:1/2048:0.5-1/2048)*80,20*log10(abs(fh8)))
grid on
axis([-40 40 -90 10])

phs=unwrap(angle(fh8));
phs=phs-phs(1025);
dphsx=conv(phs,[-1 0 1]*1024/80);
dphs=dphsx(2:2049);
dphs(1:890)=zeros(1,890);
dphs(1159:2048)=zeros(1,890);

subplot(3,1,3)
plot((-0.5:1/2048:0.5-1/2048)*80,dphs)
grid on
axis([-6 6 0 4])

figure(7)
subplot(3,1,1)
plot(0:2000,8*real(x8a(1:2001)))
hold on
plot((0:16:2001)+8,8*real(x8b(6:2:256)),'ro')
hold off
grid on
axis([0 2000 -1.5 1.5])

subplot(3,1,2)
plot((-0.5:1/4096:0.5-1/4096)*80,fftshift(20*log10(abs(fft(x8a(201:4296).*ww)))))
hold on
plot((-0.5:1/2048:0.5-1/2048)*80,fftshift(20*log10(abs(fft(h8,2048)))),'r')
hold off
grid on
axis([-40 40 -90 10])

subplot(3,1,3)
plot((-0.5:1/2048:0.5-1/2048)*10,fftshift(20*log10(abs(fft(x8b(21:2068).*w8)))))
grid on
axis([-5 5 -90 10])


%%%%%%%%%%%%%% option 4 8-path recursive Linear phase filter

reg0=zeros(1,3);
reg1=zeros(3,2);
reg2=zeros(3,2);
reg3=zeros(3,2);
reg4=zeros(3,2);
reg5=zeros(3,2);
reg6=zeros(3,2);
reg7=zeros(3,2);
v2=zeros(1,8)';
vv4=zeros(8,2125); 
m=1;
%
%x2=exp(j*2*pi*(1:17000)*1.05/8);

for n=1:8:17000-7
    xx=flipud(x2(n:n+7).');
    
    sm1a= (xx(2)-reg1(2,1))*coef3(1,1)+reg1(1,1);
    sm1b=(sm1a-reg1(3,2))*coef3(1,3)+(reg1(2,1)-reg1(3,1))*coef3(1,2)+reg1(2,2);
    
    sm2a= (xx(3)-reg2(2,1))*coef3(2,1)+reg2(1,1);
    sm2b=(sm2a-reg2(3,2))*coef3(2,3)+(reg2(2,1)-reg2(3,1))*coef3(2,2)+reg2(2,2);
    
    sm3a= (xx(4)-reg3(2,1))*coef3(3,1)+reg3(1,1);
    sm3b=(sm3a-reg3(3,2))*coef3(3,3)+(reg3(2,1)-reg3(3,1))*coef3(3,2)+reg3(2,2);
    
    sm4a= (xx(5)-reg4(2,1))*coef3(4,1)+reg4(1,1);
    sm4b=(sm4a-reg4(3,2))*coef3(4,3)+(reg4(2,1)-reg4(3,1))*coef3(4,2)+reg4(2,2);
    
    sm5a= (xx(6)-reg5(2,1))*coef3(5,1)+reg5(1,1);
    sm5b=(sm5a-reg5(3,2))*coef3(5,3)+(reg5(2,1)-reg5(3,1))*coef3(5,2)+reg5(2,2);
    
    sm6a= (xx(7)-reg6(2,1))*coef3(6,1)+reg6(1,1);
    sm6b=(sm6a-reg6(3,2))*coef3(6,3)+(reg6(2,1)-reg6(3,1))*coef3(6,2)+reg6(2,2);
    
    sm7a= (xx(8)-reg7(2,1))*coef3(7,1)+reg7(1,1);
    sm7b=(sm7a-reg7(3,2))*coef3(7,3)+(reg7(2,1)-reg7(3,1))*coef3(7,2)+reg7(2,2);
    
    v2=[reg0(3); sm1b; sm2b; sm3b; sm4b; sm5b; sm6b; sm7b];
    vv4(:,m)=fftshift(2*ifft(v2));
    m=m+1;
    
    reg0=[xx(1) reg0(1:2)];
    
    reg1(1,:)=[xx(2) reg1(1,1)];
    reg1(2,:)=[sm1a  reg1(2,1)];
    reg1(3,:)=[sm1b  reg1(3,1)];
    
    reg2(1,:)=[xx(3) reg2(1,1)];
    reg2(2,:)=[sm2a  reg2(2,1)];
    reg2(3,:)=[sm2b  reg2(3,1)];
    
    reg3(1,:)=[xx(4) reg3(1,1)];
    reg3(2,:)=[sm3a  reg3(2,1)];
    reg3(3,:)=[sm3b  reg3(3,1)];
    
    reg4(1,:)=[xx(5) reg4(1,1)];
    reg4(2,:)=[sm4a  reg4(2,1)];
    reg4(3,:)=[sm4b  reg4(3,1)];
    
    reg5(1,:)=[xx(6) reg5(1,1)];
    reg5(2,:)=[sm5a  reg5(2,1)];
    reg5(3,:)=[sm5b  reg5(3,1)];
    
    reg6(1,:)=[xx(7) reg6(1,1)];
    reg6(2,:)=[sm6a  reg6(2,1)];
    reg6(3,:)=[sm6b  reg6(3,1)];
    
    reg7(1,:)=[xx(8) reg7(1,1)];
    reg7(2,:)=[sm7a  reg7(2,1)];
    reg7(3,:)=[sm7b  reg7(3,1)];
    
end

figure(9)
for k=1:8
    subplot(2,4,k)
    plot(real(vv4(k,1:100)))
    hold on
     plot(imag(vv4(k,1:100)),'r')
    hold off
    grid on
    axis([0 100 -1.0 +1.05])
end


    
% 10 path delay line
gg=remez(239,[0 4 8 128]/128,[1 1 0 0]);
dd=20*reshape(gg,20,12);
% 
%y9=conv(y8,dd(18,:));

figure(10)
y8=vv4(5,:);
subplot(2,1,1)
plot(0:200,4*imag(y8(1:201)))
hold on
plot((0:2:500)+0,4*imag(y8((1:2:501)+0)),'ro')
hold off
grid on
axis([0 200 -1.5 1.5])

subplot(2,1,2)
plot((-0.5:1/2048:0.5-1/2048)*10,fftshift(20*log10(abs(fft(y8(1:2048).*w8)))))
grid on
axis([-5 5 -90 10])

% k
% pause(0.5)
% end