%Pramuditha A.A.H. (200476P)
%EN2063 - Signals and Systems
%Digital Filter Design - Semester 03

clc;
close all;

%Specifying design parameters.
Ap = 0.14;          %maximum passband ripple in dB.
As = 57;            %minimum stopband attenuation in dB.
wp1 = 1000;         %lower passband edge.
wp2 = 1500;         %upper passband edge.
wa1 = 700;          %lower stopband edge.
wa2 = 1700;         %upper stopband edge.
ws = 4200;          %sampling frequency.
T = 2*pi/ws;        %Sampling period.

%Selcting the transition bandwidth and cutoff frequencies.
tran_bw = min(wp1-wa1, wa2-wp2);
wc1 = wp1 - tran_bw/2;
wc2 = wp2 + tran_bw/2;

%Defining delta value using Ap and As.
delta_p = ((10^(0.05*Ap))-1)/((10^(0.05*Ap))+1);
delta_s = 10^(-0.05*As);
delta = min(delta_p, delta_s);

%Calculating the actual stopband attenuation of the filter.
Atten = -20*log10(delta);

%Designing kaiser window function. Determining alpha parameter.
if Atten<=21
   alpha = 0;
elseif Atten<=50
   alpha = (0.5842*((Atten-21)^0.4)) + (0.07886*(Atten-21));
else
   alpha = 0.1102*(Atten - 8.7);
end

%Selecting parameter D based on actual attenuation.
if Atten<=21
   D = 0.9222;
else
   D = (Atten - 7.95)/14.36;
end

%Selecting the proper odd value for N which satisfies the inequality
%N>=(ws*D/tran_bw)+1.
N = ceil((ws*D/tran_bw)+1);
if mod(N,2) == 0
   N=N+1;
end

%Forming Kaiser window wk[n].
Wk_n = zeros(N,1);
for n = -(N-1)/2:(N-1)/2
   beta = alpha * (1 - (2*n/(N-1))^2)^0.5;
   numerator = myBessel(beta);
   denominator = myBessel(alpha);
   Wk_n(n+(N-1)/2+1) = numerator/denominator;
end

%Plotting Kaiser window.
stem(Wk_n, 'filled');
title('Kaiser window');
xlabel('n');
ylabel('Wk_n');
grid on;

%Generating h(nT).
h_n = zeros(N,1);
h_n(38) = (2/ws)*(wc2-wc1);
for n = -(N-1)/2:(N-1)/2
   if n==0
      h_n(n+(N-1)/2+1) = (2/ws)*(wc2-wc1);
   else
      h_n(n+(N-1)/2+1) = (1/(n*pi)) * (sin(wc2*n*T) - sin(wc1*n*T));
   end
end
figure;

%Plotting h(nT) function.
stem(h_n, 'filled');
title('h(nT) function');
xlabel('n');
ylabel('h(nT)');
grid on;

%Getting the final digital filter.
filter = h_n.*Wk_n;

%Plotting the impulse response of the filter.
figure
stem(filter, 'filled');
title('Impulse Response of the filter');
xlabel('n + (N-1)/2');
ylabel('h[n]');
grid on;

%Plotting the magnitude response of the passband for a larger x-range.
[Amplitude, H] = freqz(filter);
analog_fr = H*ws/(2*pi);
amp_db = 20*log10(abs(Amplitude));
figure
plot(analog_fr, amp_db);
title('Magnitude Response of the filter in Passband');
xlabel('w in rad/s');
ylabel('Magnitude of H(w)');
grid on;

%Plotting the magnitude response of the passband.
[Amplitude, H] = freqz(filter);
analog_fr = H*ws/(2*pi);
amp_db = 20*log10(abs(Amplitude));
figure
plot(analog_fr, amp_db);
axis([wp1 wp2 -0.05 0.05]);
title('Magnitude Response of the filter in Passband');
xlabel('w in rad/s');
ylabel('Magnitude of H(w)');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Introducing Bessel Function to obtain zeroth order for modified bessel
%function of first kind.
function [ result ] = myBessel( x )
    k=1;
    result=0;
    term = 10;
    while (term>10^(-6))
        term = (((x/2)^k)/(factorial(k)))^2;
        result = result + term;
        k = k+1;
    end
    
    result = result+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



