%Pramuditha A.A.H. (200476P)
%EN2063 - Signals and Systems
%Digital Filter Design - Semester 03

clc;
close all;
clear;

%Specifying design parameters.
Ap = 0.14;      %Maximum passband ripple.
As = 57;        %Minimum stopband attenuation.
wp1 = 1000;     %Lower passband edge.
wp2 = 1500;     %Upper passband edge.
ws1 = 700;      %Lower stopband edge.
ws2 = 1700;     %Upper stopband edge.
fs = 4600;      %Sampling frequency. For this particular filter I had to use 4800 rad/s as my sampling rate. 

T = 2*pi/fs;    %Sampling period.

%Using prewarping technique to determine the design parameters for the
%analog filter.
Wp1 = 2*(1/T)*tan(wp1*T/2); 
Wp2 = 2*(1/T)*tan(wp2*T/2); 
Ws1 = 2*(1/T)*tan(ws1*T/2); 
Ws2 = 2*(1/T)*tan(ws2*T/2); 

%Returning filter's order and stopband cutoff frequency
[n,Ws] = cheb2ord([Wp1 Wp2]/fs,[Ws1 Ws2]/fs, Ap, As); 

disp('The order of the filter:');
disp(n)

%Returning poles, zeros and gain of nth order chebyshev type 2 filter.
[z, p, k] = cheb2ap(n, As); 

%State space representation.
[A, B, C, D] = zp2ss(z, p, k);

%Converting continuous-time state-space filter prototype to a bandpass
%filter.
[At, Bt, Ct, Dt] = lp2bp(A, B, C, D, sqrt(Wp1*Wp2), Wp2-Wp1); 

%Converting state space representation into a transfer function.
[p, q] = ss2tf(At, Bt, Ct, Dt);

%Getting the Laplase doamin frequency response.
%[h, w] = freqs(p, q, 2048);

%Implementing bilinear transformation to convert it to Z-doamin.
[Ad, Bd, Cd, Dd] = bilinear(At, Bt, Ct, Dt, (1/T));

%Getting frequency response of the digital filter in z-domain.
[hd, fd] = freqz(ss2sos(Ad, Bd, Cd, Dd), 2048, (1/T));

x = fd*4*pi^2/4800;

%Plotting the magnitude of the IIR bandpass filter in [0,pi] range.
subplot(2, 1, 1)
plot((x), 20*log10(abs(hd)))
xlim([0, pi])
title('Magnitude of IIR bandpass filter. ([0, pi] range)')
xlabel('Frequency in rad/s')
ylabel('Magnitude in dB')
grid

%Plotting the magnitude of the IIR bandpass filter in [wp1,wp2] range.
subplot(2, 1, 2)
plot(x, 20*log10(abs(hd)))
xlim([wp1,wp2]*T)
title('Magnitude of IIR bandpass filter. ([wp1, wp2] range)')
xlabel('Frequency in rad/s')
ylabel('Magnitude in dB')
grid


