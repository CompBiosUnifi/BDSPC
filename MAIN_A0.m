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

%Define a loop to scan its components and reconstruct the signal.

%Turn this operation into a function and test calling it on another column.

%Visualize multiple signals, taking into account their associated sampling frequencies.

%Provide the possibility of performing basic operations on the signal(s).

