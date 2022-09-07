% Noisy speech Audio
nBits = 16 ; 
nChannels = 2;

load noisy_speech.mp3;
[y, Fs]= audioread('noisy_speech.mp3')
nfft=length(y)

%frequency domain

xlim([100 250.9])
ylim([-0.0016 0.0044])
f=linspace(0,Fs,nfft);
y=abs(fft(y,nfft)/nfft);
freqz(nfft)

figure(1)
subplot(3,3,1)
plot(f(1:nfft/2),y(1:nfft/2))
title('frequency domain')

figure(1)
subplot(3,3,4)
plot(f(1:nfft/2),y(1:nfft/2))
title('frequency domain after filtering')