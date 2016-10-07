%%High FRET histogram (ALL), quick script for assembling all the
%%FRET time traces of intermediate FRET analysed molecules - RC 15/07/12
clear; close all; clc

loadnames = uigetfile({'*.mat'}, 'Select analysed files','MultiSelect', 'on');
E_data_all = [];

for i = 1:length(loadnames)
    filename = loadnames{i};
    fprintf(['\nNow analysing: \n' filename '\n'])
    load(filename)
    E_data_all = [E_data_all; E_b4_bleach];
    clear E_b4_bleach;
end

save('FRET_ALL_120820.mat');

rect = [500, 500, 350, 350]
figure('OuterPosition',rect)
hist(E_data_all,0:0.02:1);xlim([0 1])