clear all
close all

%notch filter design
Fs=16e3;
Ts=1/Fs;

M=128;%order
N=M+1;%no of filter coefficients
fo=1e03*Ts;%notch freq1
f1=2e03*Ts;%notch freq2
f2=4e03*Ts;%notch freq3

wo=2*pi*fo;
w1 =2*pi*f1;
w2 =2*pi*f2;

%notch width
alfa =10*2*pi*Ts;

%FInd the P matrix
P=zeros(N);
q=zeros(N,1);

%part1
dw=pi/400;%freq domain sampling
w=[0:dw:(wo-alfa/2)]';
nM=[0:N-1]';
U=cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*trapz(U,2)*dw;

%part2
dw=pi/1000;

W=500;%Weight of Notch part
epsi=0.001;
w=[(wo-alfa/2):dw:(wo+alfa/2)]';
nM=[0:N-1]';
U=cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+W*trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*W*epsi*trapz(U,2)*dw;

%part3
dw=pi/400;
w=[(wo+alfa/2):dw:(w1-alfa/2)]';
nM=[0:N-1]';
U = cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+trapz(U(n1,:).*U(n2,:))*dw;
end

end
q=q-2*trapz(U,2)*dw;

%part4
dw=pi/1000;
W=500;%Weight of Notch part
epsi=0.001;
w=[(w1-alfa/2):dw:(w1+alfa/2)]';
nM=[0:N-1]';
U=cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+W*trapz(U(n1,:).*U(n2,:))*dw;

end
end
q=q-2*W*epsi*trapz(U,2)*dw;

%part5
dw=pi/400;
w=[(w1+alfa/2):dw:(w2-alfa/2)]';
nM=[0:N-1]';
U=cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+W*trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*trapz(U,2)*dw;

%part6
dw=pi/1000;
W=500;%Weight of Notch part
epsi=0.001;
w=[(w2-alfa/2):dw:(w2+alfa/2)]';
nM=[0:N-1]';
U=cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+W*trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*W*epsi*trapz(U,2)*dw;

%part7
dw=pi/400;
w=[(w2+alfa/2):dw:pi]';

nM=[0:N-1]';
U=cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+W*trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*trapz(U,2)*dw;

figure(1);
plot(q)

%solve for minimization
a=-P\q;
for k=1:M/2-1
h(M/2-k)=a(k+1)/2;
h(M/2+k)=a(k+1)/2;
end
h(M/2)=a(1);
plot(h);
figure(2)
stem(h),shg
figure(3)
F=linspace(0,Fs/2,2000);
H=freqz(h,1,F,Fs);
plot(F,abs(H)),grid on

This is only for 2 notchs

clc
clear all
close all
%comb filter design
Fs=10000;

Ts=1/Fs;

M=128;%order
N=M/2+1;%no of filter coefficients

%comb width
alfa=10*2*pi*Ts;

fo=1000*Ts;%notch freq
wo=2*pi*fo;

f1 = 2000*Ts;
w1 = 2*pi*f1;

%FInd the P matrix
P=zeros(N);
q=zeros(N,1);

%part1
dw=pi/300;%freq domain sampling
w=[0:dw:(wo-alfa/2)]';
nM=[0:N-1]';
U=cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*trapz(U,2)*dw;

%part2

dw=pi/500;
W=10000;%Weight of Notch part
epsi=0.001;
w=[(wo-alfa/2):dw:(wo+alfa/2)]';
nM=[0:N-1]';
U=cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+W*trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*W*epsi*trapz(U,2)*dw;

%part3
dw=pi/300;
w=[(wo+alfa/2):dw:(w1-alfa/2)]';
nM=[0:N-1]';
U = cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*trapz(U,2)*dw;
%figure(1);
%plot(q)

%part 4

dw=pi/500;
w=[(w1-alfa/2):dw:(w1+alfa/2)]';
nM=[0:N-1]';
U = cos(nM*w');
for n1=1:N
for n2=1:N

P(n1,n2)=P(n1,n2)+W*trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*W*epsi*trapz(U,2)*dw;

%part 5
dw=pi/300;
w=[(w1+alfa/2):dw:pi]';
nM=[0:N-1]';
U = cos(nM*w');
for n1=1:N
for n2=1:N
P(n1,n2)=P(n1,n2)+trapz(U(n1,:).*U(n2,:))*dw;
end
end
q=q-2*trapz(U,2)*dw;

%solve for minimization
a=-P\(q/2);
for k=1:M/2-1
h(M/2-k)=a(k+1)/2;
h(M/2+k)=a(k+1)/2;
end
h(M/2)=a(1);

figure(2)
stem(h),shg
title('Impulse Response')
xlabel('N')

ylabel('h[n]')
figure(3)
subplot(211)
F=linspace(0,Fs/2,2000);
H=freqz(h,1,F,Fs);
plot(F,abs(H));
xlabel('Freqeuncy')
ylabel('Magnitude')
title('Freq Response')
grid on
subplot(212)
plot (F, phase(H));
title('Phase Response')