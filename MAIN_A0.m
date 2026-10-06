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
L5 – ECG/HRV Signal Analysis – Part B

%subwindowing the entire signal into windows of N samples with a defined overlap
% Define window size and overlap for subwindowing
fs=256;
windowSize = 60*fs; % Number of samples in each window
overlap = windowSize/2;    % Number of overlapping samples

% Create overlapping windows for the signal
numWindows = floor((length(all_sig) - windowSize) / overlap) + 1; %<<-- valid only for 50% overlap
windows_ecg = cell(2, numWindows);

%extract the windows from the signal and store them in a cell array
for j = 1:numWindows
    startIdx = (j-1) * overlap + 1;
    windows_ecg{1, j}= all_sig(startIdx:startIdx + windowSize - 1);
    windows_ecg{2, j}= [startIdx startIdx + windowSize - 1];
end

%feature extraction for each segmented window
for j=1:numWindows
    curr_wind=windows_ecg{1,j};
    %PT and RR
    [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(curr_wind,fs,0);
    RR=diff(qrs_i_raw/fs); %obtain the tachogram signal in [s]
    RR_round = round(RR*fs)/fs;
    
    %Discussion point 1: do all windows provide good R-peak detection?
    %Discussion point 2: how can possible anomalies be identified? Before
    %or after feature extraction?
    
    %example of feature extraction (rough implementation!) -
    %Discussion point 3: is this the most efficient method? I need to
    %keep track of the feature names somewhere
    feat_all(1,j)=HRV.HR(RR,0);
    feat_all(2,j)=HRV.RMSSD(RR*1000,0);
    [pLF,pHF,LFHFratio,VLF,LF,HF,~,~,~,TP] = HRV.fft_val_fun(RR,fs,'linear');
    feat_all(3,j)=LFHFratio;

end

%HR 
% hr = HRV.HR(RR,0); %0 is the entire recording
% %RMSSD
% rmssd = HRV.RMSSD(RR*1000,0); %conversion to ms 
% %LF - HF
% [pLF,pHF,LFHFratio,VLF,LF,HF,~,~,~,TP] = HRV.fft_val_fun(RR,fs,'linear');

%turn this operation (feature extraction) into a function with
%user-defined parameters.

%Produce a plot summarizing the trend of a series of these metrics as a function of the window.
%Plot the metrics as a function of the window

figure;
% windowTime = (0:numWindows-1) * (windowSize/2) / fs;

subplot(3,1,1)
plot(1:numWindows, feat_all(1,:), '-o', 'LineWidth', 1.2);
xlabel('Window');
ylabel('HR [bpm]');
title('Heart Rate');
grid on;

subplot(3,1,2)
plot(1:numWindows, feat_all(2,:), '-o', 'LineWidth', 1.2);
xlabel('Window');
ylabel('RMSSD [ms]');
title('RMSSD');
grid on;

subplot(3,1,3)
plot(1:numWindows, feat_all(3,:), '-o', 'LineWidth', 1.2);
xlabel('Window');
ylabel('LF/HF');
title('LF/HF ratio');
grid on;

%Consider the window size and its consistency in HRV metric extraction.
%Discussion point: window size, extractable features, overlap.

%%% >--------------------------------------------------------------------------------<
L6 – Characterization of HRV States and Events – Part C

% -- > code for labeling periods from the .txt file
% -- > sample-wise assignment of labels from the .txt file

% e.g. given the start and end sample, if it falls within this range <--
% this information also needs to be returned by the segmentation function.

sleep_label = import_sleep_info('SN001_sleepscoring.txt');

%how to handle conflicts, i.e., transitions?

%e.g. 1 grouping
dict={"Sleep stage W", "Sleep stage N1", "Sleep stage N2", "Sleep stage N3", ...
    "Sleep stage R","Extra_label";...
    0 1 2 3 4 9};
annotat_sleep=sleep_label.Annotation;
onset_sleep=sleep_label.RecordingOnset;

for i=1:size(windows_ecg,2)
    current_samples=windows_ecg{2,i}./fs; %<<-- in seconds
    %find the corresponding label 
    idx_good=onset_sleep >= current_samples(1,1) & onset_sleep <= current_samples(1,2);
    %candidate/s?
    cand_label=annotat_sleep(idx_good);
    if size(cand_label,1) == 1
        the_label=strcmp([dict{1,:}],cand_label);
        if ~isempty(the_label)
            wind_label_ecg(1,i)=dict{2,the_label};
        else
             wind_label_ecg(1,i)=dict{2,end};
        end

    else % a "decision rule" is needed — "majority"?
        %define the majority rule on cand_label to find the winner
        %here
        % Apply majority voting to determine the final label
        for j=1:size(dict,2)
            count_label(j)=sum(contains(cand_label,dict{1,j}));
        end
              %what if there are ties? --- handle this exception!
        %e.g. 1A (take only the last observed transition)
        [~, maxIdx] = max(count_label); % Find the index of the maximum count
        num_max=sum(count_label == max(count_label));

        if num_max > 1
            the_label=strcmp([dict{1,:}],cand_label(end,1)); %<<< an exception may occur here!
            %i.e., the last label may not correspond to either of
            %the two majority cases!
              wind_label_ecg(1,i)=dict{2,the_label};
        else
             wind_label_ecg(1,i) = dict{2,maxIdx};
        end
  
    end
end

histogram(wind_label_ecg(end,:));
xlabel("wind_label_ecg(end,:)");
title("wind_label_ecg(end,:)");
legend("show");

%Discussion point: limitations and fragility of the labeling rules —
%they have an impact on the system's performance

%Discussion point: how could this information be extended from an
%offline scenario to an online scenario?

%%% >--------------------------------------------------------------------------------<

