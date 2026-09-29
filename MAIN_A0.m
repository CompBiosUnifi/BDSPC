%% L0 - Open and discuss the .edf file
clear all
close all
clc

%open the edf file SN001
edfFile = 'SN001.edf';
edfInfo = edfinfo(edfFile);

% Read the data from the EDF file
edfData = edfread(edfFile);

%Punto di discussione, riesco a risalire e leggere tutte le informazioni necessarie dei file salvati?
%Can I access and read all the necessary information from the saved files?

%%% >--------------------------------------------------------------------------------<
%% L1 - Open and discuss the content, and visualize the entire dataset.
%Different ways to accomplish the same thing.

%Open and read the .edf file used in L0.
edfFile = 'SN001.edf';
edfInfo = edfinfo(edfFile);

% Read the data from the EDF file
edfData = edfread(edfFile);

%Define a loop to scan its components and reconstruct the signal.
idx_cl=8;
all_sig=[];
for i=1:5000 %size(edfData,1) 
    all_sig=[all_sig; edfData.(idx_cl){i,:}];
end

figure(1)
plot(all_sig)

% TO DO AT HOME
%Turn this operation into a function and test calling it on another column.
%Visualize multiple signals, taking into account their associated sampling frequencies.
%Provide the possibility of performing basic operations on the signal(s).

%%% >--------------------------------------------------------------------------------<
%% L2 - DFT application

%Spectrum visualization on segments of the extracted signal.
% FFT on a portion of the signal
x_signal=all_sig(15900:15900+256*3+1);
Fs = 256;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(x_signal);             % Length of signal
% t = (0:L-1)*T;        % Time vector
Y = fft(x_signal);
P2 = abs(Y); % non-normalized coefficients for aperiodic signals (/L)  % non-normalized coefficients for aperiodic signals (the DFT summation should be considered)
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1); %(correction required when excluding the negative-frequency components)

figure(2)
f = Fs*(0:(L/2))/L;
plot(f,P1) 
% hold on
title('Single-Sided Amplitude Spectrum of xsignal')
xlabel('f (Hz)')
ylabel('Amp xsignal')
axis tight

%Evaluation of increasing the number of DFT points.
n=length(x_signal)*10; 
% n = 2^nextpow2(length(x_signal)); % for efficient FFT rather than DFT.

Fs = 256;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
L = length(x_signal);             % Length of signal
% t = (0:L-1)*T;        % Time vector
Y = fft(x_signal,n);

P2 = abs(Y); % non-normalized coefficients for aperiodic signals
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(n/2))/n;
figure(6)
plot(f,P1) 
% hold on
title('Single-Sided Amplitude Spectrum of xsignal - lower number of FFT points')
xlabel('f (Hz)')
ylabel('Amp xsignal')
axis tight

%Discussion points: 
%Are these representations consistent with the a priori information known about the signals?
%Is the representation window correct?

%%% >--------------------------------------------------------------------------------<
%%L3 - Extraction of the spectrum from the entire signal, its segments, and visualization
% Identify the properties of the sampled signal
fs = 256;            % Sampling frequency                    
T = 1/fs;             % Sampling period       
N = length(all_sig);  % Length of signal
t = (0:N-1)*T;        % Time vector

[pxx_welch,f] = pwelch(all_sig,N/100,[],512,fs,'psd');
pwelch(all_sig,N/2,100,512,fs,'psd');

% Hypothesis that different phases carry non-stationary frequency information

% Isolation of these phases
N1_samp=all_sig(240*fs:240*fs+10*fs-1);
N=length(N1_samp);
figure(20)
pwelch(N1_samp,N/10,[],512,fs,'psd');
hold on

N3_samp=all_sig(3360*fs:3360*fs+10*fs-1);
N=length(N3_samp);
pwelch(N3_samp,N/10,[],512,fs,'psd');

R_samp=all_sig(4610*fs:4610*fs+10*fs-1);
N=length(R_samp);
pwelch(R_samp,N/10,[],512,fs,'psd');

hold off

% Discussion point: does the window remain consistent for these processes?

%%% >--------------------------------------------------------------------------------<
%% L4 - ECG/HRV Signal Analysis Part A

% display the extracted ECG signal
figure(2)
plot(all_sig)
fs=256;

%all_sig(1)=all_sig(1)+10000000;

% implement R-peak detection % Pan-Tompkins used by user-function from mathworks repository
[qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(all_sig,fs,1);
RR=diff(qrs_i_raw/fs); % obtain the signal for the tachogram in [s]
RR_round = round(RR*fs)/fs;
Ann = cumsum(RR_round);
figure(12)
plot(Ann,RR_round)
ylabel('RR (s)')
xlabel('Time (s)')
title('Tachogram')

% turn this operation (detection and visualization) into a function
%TO DO AT HOME

% Discussion point: tachogram and its alterations in the observed signal

%%% >--------------------------------------------------------------------------------<


