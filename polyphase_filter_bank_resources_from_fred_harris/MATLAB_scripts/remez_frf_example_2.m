% Remez filter design example
% Design a low pass FIR filter with the following characteristics
% filter length is a multiple of 30
% sample rate 4.5 MHz
% passband edge 50 kHz
% passband ripple 0.1 dB
% stopband edge 100 kHz
% rejection 60 dB
% the weights are adjusted to optimize the inband vs. stopband performance
h1=remez(280-1,[0, 1, 2, 50]/50, [1 1 0 0], [1 10]);

figure(1)
plot((-0.5:1/4096:0.5-1/4096),20*log10(abs(fftshift(fft(h1,4096)))))
axis([-0.5 0.5 -100 5])
title('Equal Ripple Remez Filter Design Example')
grid on
zoom on
% now give the filter a 1/f rolloff in the stopband
% myfrf is a modified version of the remezfrf function to provide the 1/f rolloff
%   in myfrf, removed linse start with "%---"
%   and new lines end with "%+++ ..."
%   to make the changes visible
% the weights are adjusted to optimize the inband vs. stopband performance
h2=remez(280-1,[0, 1, 2, 50]/50, {'myfrf', [1 1 0 0]}, [1 10]);
figure(2)
plot((-0.5:1/4096:0.5-1/4096),20*log10(abs(fftshift(fft(h2,4096)))))
axis([-0.5 0.5 -100 5])
title('1/f Stopband Rolloff Remez Filter Design Example')
grid on
zoom on
