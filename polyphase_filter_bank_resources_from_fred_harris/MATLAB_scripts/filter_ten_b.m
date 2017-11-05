%filter_ten_a

% Animated spectra and time response of 10-stage polyphase filter bank.
% Can move single tone through filter bank via a slider control or by a scheduled sweep.
% can view time series in each channel or can view spectrum in each channel.
% program starts by depressing one of 4 GUI push buttons. 
% This file calls filter_ten_b_call
% scrip written by fred harris of SDSU
% bandwidth greater than specral spacing

clear all
close all

freq=0.1;
 xx= exp(j*2*pi*(0:2499)*freq);

figure(1)
subplot(6,2,1)
d_in=zeros(1,100);
plot(0:99,real(d_in));
grid
axis([0 100 -1.1 1.1]);

subplot(6,2,3)
d_1=zeros(1,100);
plot(0:0.1:9.9,d_1);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,4)
d_2=zeros(1,100);
plot(0:0.1:9.9,d_2);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,5)
d_3=zeros(1,100);
plot(0:0.1:9.9,d_3);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,6)
d_4=zeros(1,100);
plot(0:0.1:9.9,d_4);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,7)
d_5=zeros(1,100);
plot(0:0.1:9.9,d_5);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,8)
d_6=zeros(1,100);
plot(0:0.1:9.9,d_6);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,9)
d_7=zeros(1,100);
plot(0:0.1:9.9,d_7);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,10)
d_8=zeros(1,100);
plot(0:0.1:9.9,d_8);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,11)
d_9=zeros(1,100);
plot(0:0.1:9.9,d_9);
grid
axis([0 10 -1.1 1.1]);

subplot(6,2,12)
d_10=zeros(1,100);
plot(0:0.1:9.9,d_10);
grid
axis([0 10 -1.1 1.1]);

freq=0.1;
flag1=1;    
% flag1=1 slider control, flag1=0 scheduled sweep
flag2=1;    
% flag2=1 time display, flag2=0 frequency display

slider_1=uicontrol('style','slider','units','normalized','pos',[0.65 0.90 0.19 0.028],...
         'min',-0.4,'max',+0.4,'value',freq,...
         'callback',['freq=0.01*round(100*get(slider_1,''value''));',...
             'set(slider_1_cur,''string'',num2str(freq)),',...
             'set(gca,''view'',[0 0.5]),',...
             'filter_ten_b_call(freq,flag1,flag2)']);

slider_1_min=uicontrol('style','text','units','normalized','pos',[0.62 0.90 0.025 0.025],...
    'string', num2str(get(slider_1,'min')));

slider_1_max=uicontrol('style','text','units','normalized','pos',[0.84 0.90 0.025 0.025],...
    'string',num2str(get(slider_1,'max')));

slider_1_cur=uicontrol('style','text','units','normalized','pos',[0.78 0.87 0.055 0.025],...
    'string',num2str(0.01*round(100*get(slider_1,'value'))));

slider1_title=uicontrol('style','text','units','normalized','pos',[0.690 0.87 0.10 0.025],...
    'string','Center Frequency');


h30=uicontrol('style','pushbutton','string','SLIDER','units','normalized',...
    'position',[0.6 0.830 0.055 0.030],'callback',['flag1=1;','filter_ten_b_call(freq,flag1,flag2)']);
h40=uicontrol('style','pushbutton','string','SWEEP','units','normalized',...
    'position',[0.67 0.830 0.055 0.030],'callback',['flag1=0;','filter_ten_b_call(freq,flag1,flag2)']);
h50=uicontrol('style','pushbutton','string','TIME','units','normalized',...
    'position',[0.74 0.830 0.055 0.030],'callback',['flag2=1;','filter_ten_b_call(freq,flag1,flag2)']);
h60=uicontrol('style','pushbutton','string','FREQ','units','normalized',...
    'position',[0.81 0.830 0.055 0.030],'callback',['flag2=0;','filter_ten_b_call(freq,flag1,flag2)']);

% gg=get(gca);
% set(gca,'gridlinestyle','-')