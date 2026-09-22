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

%Evaluation of increasing the number of DFT points.

%Discussion points: 
%Are these representations consistent with the a priori information known about the signals?
%Is the representation window correct?


%%% >--------------------------------------------------------------------------------<
