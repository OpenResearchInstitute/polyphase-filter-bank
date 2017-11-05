%Polyphase Channelizer demo

clear all;
close all;

mhz = 1e6;

%Create reference FSK signal at Fs/4
Fs = 40e6;
Rb = 500e3;
xb = floor(2*rand(1,400));
sps = Fs/Rb;
ref_FSK = fskmod(xb,2,Rb,sps,Fs);
n = 0:length(ref_FSK)-1;
ref_FSK = real(ref_FSK.*exp((-j*2*pi/4).*n))*2^11;

%Create baseband noise
noise=randn(1,length(ref_FSK));

%Create interferer at 19 Mhz
%intf = cos(2*pi*19e6/Fs.*n)*2^11
intf = (cos(2*pi*19e6/Fs.*n)+cos(2*pi*18e6/Fs.*n)+cos(2*pi*17e6/Fs.*n))*2^11;
intf=4*intf;

%Design a quick and dirty filter to make shaped noise
Fc = 2e6;
b = fir1(250,Fc/(Fs/2));
filtnoise = filter(b,1,noise)*2^11;

ref_ADC = ref_FSK+filtnoise+intf;

figure(1);
Nfft = length(ref_ADC);
H = abs(fft(ref_ADC,Nfft));
H = 20*log10(H./max(H));
plot([-.5:1/Nfft:.5-1/Nfft]*Fs/mhz,fftshift(H));
title('Four Channel Input Spectrum');
xlabel('Freq (MHz)');
ylabel('Magnitude (dB)');
axis([-Fs/(2*mhz) Fs/(2*mhz) -90 5]);
grid on

pause 

% define channel filter specifications
polyChan.N = 160;
polyChan.M = 4;
polyChan.Fs = Fs;
polyChan.Fpass = 3e6;
polyChan.Fstop = 5e6;
polyChan.CL = 16;
polyChan.FL = 15;
polyChan.Atten_dB = 80;
polyChan.N_est = (polyChan.Fs/(polyChan.Fstop-polyChan.Fpass))*polyChan.Atten_dB/22; 

% design filter
polyChan.f = [0 polyChan.Fpass polyChan.Fstop polyChan.Fs/2]/(polyChan.Fs/2);       %frequency vector
polyChan.a = [1 1 0 0];                             %amplitude vector
polyChan.w = [1 2000];                              %weight vector
polyChan.b=remez(polyChan.N-1,polyChan.f, {'myfrf',polyChan.a},polyChan.w);
polyChan.b = polyChan.b/max(polyChan.b);

%quantize filter
q = quantizer('fixed','round','saturate',[polyChan.CL polyChan.FL]);
polyChan.b_q = quantize(q,polyChan.b);
polyChan.b_coef = polyChan.b_q*2^polyChan.FL;

bb = zeros(4,polyChan.N);
for k=1:4
    bb(k,k:4:polyChan.N) = polyChan.b_q(k:4:polyChan.N);
end

% plot phase rocket
figure(2)
plot(0,0)
hold on
Nfft = polyChan.N;
for nn=1:4
    phs=unwrap(angle(fftshift(fft(bb(nn,:),Nfft))));
    phs_sv(nn)=phs((Nfft/2)+1);
    plot([-0.5:1/Nfft:0.5-1/Nfft]*20,phs-phs_sv(nn),'linewidth',2)
end
grid on
title('Spectral Phase Profiles across Nyquist Zones -2 to 2')

pause

%compute reference phase
Nfft = polyChan.N;
bb1=zeros(1,polyChan.N);
bb_nc=fftshift(polyChan.b_q);
bb1(1:4:polyChan.N)=bb_nc(1:4:polyChan.N);
fbb1=fftshift(abs(fft(bb1,Nfft)));
arg1=unwrap(fftshift(angle(fft(bb1,Nfft))))/(2*pi);
arg_r=arg1;

%plot channelizer filter showing phase profiles in 3-space
figure(3);
for m=1:4
    subplot(2,2,m);
    bb1=zeros(1,polyChan.N);
    bb1(m:4:polyChan.N)=bb_nc(m:4:polyChan.N);
    fbb1=fftshift(abs(fft(bb1,Nfft)));
    arg1=unwrap(fftshift(angle(fft(bb1,Nfft))))/(2*pi);
    arg1=arg1-arg_r;
    plot3(imag(fbb1.*exp(j*2*pi*arg1)), [-0.5:1/Nfft:0.5-1/Nfft]*20,real(fbb1.*exp(j*2*pi*arg1)),'r','linewidth',2);
    title(['Spectral Phase Profile for sub-filter ',num2str(m-1)]);
    xlabel('Imag');
    ylabel('Freq');    
    zlabel('Real');    
    grid on
    %axis([-1.5 1.5 -0.5 0.5  -1.5 1.5])
        axis([-1.5 1.5 -10 10  -1.5 1.5])

    view([110,20])
    line([0 0],[0 0],[-1.2 1.2],'linewidth',2); 
    %line([0 0],[-.5 .5],[0 0],'linewidth',2);
   line([0 0],[-10 10],[0 0],'linewidth',2);
   
    line([-1.2 1.2],[0 0],[0 0],'linewidth',2);
    
end

pause;

%Overlay input spectrum with Channelizer freq response - 
%artifically mixed up to 10Mhz to show what is equivalently happening.
figure(4);
plot([-.5:1/length(H):.5-1/length(H)]*Fs/mhz,fftshift(H));
title('ADC input spectrum overlayed with Channelizer response');
xlabel('Freq (MHz)');
ylabel('Magnitude (dB)');
grid on
hold on
Nfft = 1024;
H = 20*log10(abs(fft(bb(1,:)./sum(bb(1,:)),Nfft)));
plot([-.5:1/Nfft:.5-1/Nfft]*Fs/mhz,fftshift(H(1,:)),'-r','linewidth',2);
axis([-Fs/(2*mhz) Fs/(2*mhz) -90 5]);

pause 

%Run the channelizer in slow motion...

%Create polyphase partition of filter
bb = reshape(polyChan.b_q,4,polyChan.N/4);

%Partition the input data as 4 4:1 downsampled streams
xx = reshape(ref_ADC,4,length(ref_ADC)/4);

%run downsampled streams through the polyphase filter arms
yy(4,:) = filter(bb(4,:),1,xx(1,:));
yy(3,:) = filter(bb(3,:),1,xx(2,:));
yy(2,:) = filter(bb(2,:),1,xx(3,:));
yy(1,:) = filter(bb(1,:),1,xx(4,:));

%sum down each column
y = sum(yy);  

figure(5);
Nfft = length(y);
H = abs(fft(y,Nfft));
H = 20*log10(H./max(H));
plot([-.5:1/Nfft:.5-1/Nfft]*Fs/(4*mhz),fftshift(H));
title('Output Spectrum resulting from 4-Fold aliasing to baseband');
xlabel('Freq (MHz)');
ylabel('Magnitude (dB)');
axis([-Fs/(8*mhz) Fs/(8*mhz) -90 5]);
grid on

Nsteps = 50;
%%%%% Select the channel here
k = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 0:1/Nsteps:1;
rot(4,:) =  exp(j*2*pi*(3/4)*k.*n);
rot(3,:) =  exp(j*2*pi*(2/4)*k.*n);
rot(2,:) =  exp(j*2*pi*(1/4)*k.*n);
rot(1,:) =  ones(1,length(n));
 
pause 

bb1=zeros(4,polyChan.N);
for m=1:4
    arg_ref(m,:)=arg_r;
end

% Plot all subfilters overlayed as a individual figure first
figure(6);
for nn=1:4
    rot1(nn) = exp(j*2*pi*(nn-1)*m*k/4);
    bb1(nn,nn:4:polyChan.N)=bb_nc(nn:4:polyChan.N)*rot1(nn);
end

Nfft = polyChan.N;
fbb1=fftshift(abs(fft(bb1,Nfft,2)));
arg1=unwrap(fftshift(angle(fft(bb1,Nfft,2))))/(2*pi);
arg1=arg1-arg_ref;
plot3(imag(fbb1.*exp(j*2*pi*arg1)), [-0.5:1/Nfft:0.5-1/Nfft]*20,real(fbb1.*exp(j*2*pi*arg1)),'r','linewidth',1.5);
title('Overlay of Spectral responses');
xlabel('Imag');
ylabel('Freq');    
zlabel('Real');    
grid on

line([0 0],[0 0],[-1.2 1.2],'linewidth',1.5);
line([0 0],[-10 10],[0 0],'linewidth',1.5);
line([-1.2 1.2],[0 0],[0 0],'linewidth',1.5);

axis([-1.5 1.5 -10 10  -1.5 1.5])
view([110,20])

pause;


hdl=figure;
zz=1;
for m=0:1/Nsteps:1

    % Show filters as the alias separates
    for nn=1:4
        rot1(nn) = exp(j*2*pi*(nn-1)*m*k/4);
        bb1(nn,nn:4:polyChan.N)=bb_nc(nn:4:polyChan.N)*rot1(nn);
    end

    Nfft = polyChan.N;
    fbb1=fftshift(abs(fft(bb1,Nfft,2)));
    arg1=unwrap(fftshift(angle(fft(bb1,Nfft,2))))/(2*pi);
    arg1=arg1-arg_ref;
    subplot(2,2,1)
    plot3(imag(fbb1.*exp(j*2*pi*arg1)), [-0.5:1/Nfft:0.5-1/Nfft]*20,real(fbb1.*exp(j*2*pi*arg1)),'r','linewidth',1.5);
    title(['Separating the k=',num2str(k),' Alias - Iteration ',num2str(m*Nsteps)]);
    xlabel('Imag');
    ylabel('Freq');    
    zlabel('Real');    
    grid on

    line([0 0],[0 0],[-1.2 1.2],'linewidth',1.5);
    line([0 0],[-10 10],[0 0],'linewidth',1.5);
    line([-1.2 1.2],[0 0],[0 0],'linewidth',1.5);
 
    %tick marks
%     line([0 0 0 0 0 0 0 0 0],[1  1  1  2  2  2  3  3  3]/4,[0 .05 0 0 .05 0 0 .05 0]);
%     line([0 0 0 0 0 0 0 0 0],-[ 1  1  1  2  2  2  3  3  3]/4,[0 .05 0 0 .05 0 0 .05 0]);

%     text(0, 0.58, 0,'Freq','fontsize',10,'fontweight','bold')
%     text(0, -0.02, 1.4,'Real','fontsize',10,'fontweight','bold')
%     text(1.6, 0, 0,'Imag','fontsize',10,'fontweight','bold')

    axis([-1.5 1.5 -10 10  -1.5 1.5])
    view([110,20])

%     xlabel('imaginary part','rotation',15)
%     ylabel('frequency          ','rotation',-8)
%     zlabel('real part')
        
    % display phase coherent sum
    fbb1=fftshift(fft(bb1,Nfft,2));
    fbb1_sum = abs(sum(fbb1))/4;
    arg1=unwrap(angle(sum(fbb1)/4))/(2*pi);
    arg1=arg1-arg_r;
    subplot(2,2,2)
    plot3(imag(fbb1_sum.*exp(j*2*pi*arg1)), [-0.5:1/Nfft:0.5-1/Nfft]*20,real(fbb1_sum.*exp(j*2*pi*arg1)),'r','linewidth',1.5);
    xlabel('Imag');
    ylabel('Freq');    
    zlabel('Real');    
    grid on
    
    line([0 0],[0 0],[-1.2 1.2],'linewidth',1.5);
    line([0 0],[-10 10],[0 0],'linewidth',1.5);
    line([-1.2 1.2],[0 0],[0 0],'linewidth',1.5);
 
    %tick marks
%     line([0 0 0 0 0 0 0 0 0],[1  1  1  2  2  2  3  3  3]/4,[0 .05 0 0 .05 0 0 .05 0]);
%     line([0 0 0 0 0 0 0 0 0],-[ 1  1  1  2  2  2  3  3  3]/4,[0 .05 0 0 .05 0 0 .05 0]);

%     cc=get(gca,'children');
%     set(cc,'linewidth',1.5)

%     text(0, 0.58, 0,'Freq','fontsize',10,'fontweight','bold')
%     text(0, -0.02, 1.4,'Real','fontsize',10,'fontweight','bold')
%     text(1.6, 0, 0,'Imag','fontsize',10,'fontweight','bold')

    title('Phase Coherent Sum');
    axis([-1.5 1.5 -10 10  -1.5 1.5])
    view([110,20])

%     xlabel('imaginary part','rotation',15)
%     ylabel('frequency          ','rotation',-8)
%     zlabel('real part')
    
    % Show Data as the alias separates
    for mm=4:-1:1
        index = round((m*Nsteps)+1);
        yy1(mm,:) = yy(mm,:)*rot(mm,index);
    end
    subplot(2,2,4)
    y = sum(yy1);
    Nfft = length(y);
    H = abs(fft(y,Nfft));
    H = 20*log10(H./max(H));
    plot([-.5:1/Nfft:.5-1/Nfft]*Fs/(4*mhz),fftshift(H));
    axis([-Fs/(8*mhz) Fs/(8*mhz) -90 5]);
    grid on
    title('Channelizer Output Spectrum');
    subplot(2,2,3)
    compass(real(rot1),imag(rot1));
    title('Phase Rotator Progression')
    
    if( m==0 )
        pause; 
    else
        pause(.1)
    end

%UNCOMMENT the following 2 lines and the movie2avi() line to create an avi movie
%     
%      FF(zz)=getframe(hdl);
%      zz=zz+1;
%      FF(zz)=getframe(hdl);
%      zz=zz+1;
  end
% pause
%  movie2avi(FF,'paddlewheel_2');
