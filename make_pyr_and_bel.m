% This scipt can be used to obtain the free artifact downsampled signals
% from the CBD dataset batch 1 and 2

% Script designed to run the analysis of SWRs
% Load various folders with relevants functions afterward
clear
clc

% SET THE ANIMAL NUMBER HERE
% __________
animal = "209";
batch = "2";
% ----------

source_path = uigetdir('','Source Path');
results_folder = fullfile(source_path, "sleep_oscil_tools", "results", "artifact_free");
if ~exist(results_folder, 'dir')
    mkdir(results_folder)
end

% Sleep Oscilations Tools (self)
addpath(fullfile(source_path, "sleep_oscil_tools"));

% Rat N - vehicle
data_path = fullfile(source_path, "data");

%Paths
sleep_states_path = fullfile(data_path, "PCA_state_files", ...
                          "PCA_Scorer_batch_" + batch, ...
                          "Rat" + animal, "states.mat" );
pyra_path = fullfile(data_path, "HPC_all_animals", ...
                "HPC_" + animal + "_CHpyr.continuous.mat");
abov_path = fullfile(data_path, "HPC_all_animals", ...
                "HPC_" + animal + "_CHab.continuous.mat");
belo_path = fullfile(data_path, "HPC_all_animals", ...
                "HPC_" + animal + "_CHbel.continuous.mat");
shal_path = fullfile(data_path, "PFC_all_animals", ...
                "PFC_" + animal + "_CH7.continuous.mat");
deep_path = fullfile(data_path, "PFC_all_animals", ...
                "PFC_" + animal + "_CH19.continuous.mat");

% Acquisition parameters
%acq_fhz = 30e3; %acquisition freq
fn = 600;   %downsampling freq

% Sleep states across the signal
clusters = importdata(sleep_states_path);
sleep_states = zeros(1, length(clusters)*fn*10);
for i = 1:length(clusters)
    state = repmat(clusters(i),1,6000);
    pos = 6000*i-5999;
    sleep_states(pos:pos+5999) = state;
end
L = length(sleep_states);

% Load data 
pyra = importdata(pyra_path);
abov = importdata(abov_path);
belo = importdata(belo_path);
shal = importdata(shal_path);
deep = importdata(deep_path);

% Crop data to the same length
pyra = pyra(1:L);
abov = abov(1:L);
belo = belo(1:L);
shal = shal(1:L);
deep = deep(1:L);

% Sum of amplitudes to visually assess artifacts as they will appear in
% every channel and add up
TOT = abs(pyra) + abs(abov) + abs(belo) + abs(shal)+ abs(deep);

%figure
threshold = 4500;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), TOT)
yline(threshold, '-', "Threshold:" + int2str(threshold))
title("Channels Amplitudes Sum: Rat " + animal)
%%
% Artifacts
tr = 5670; %Visual threshold 
outliers = false(L,1);
index = 1;
while index<L
    if TOT(index)>=tr
        if index-300 < 1; lo_lim = 1; else; lo_lim = index-300; end
        if index+1999 > L; hi_lim = L; else; hi_lim = index+1999; end
        sz = length(lo_lim:hi_lim);
        outliers(lo_lim:hi_lim) = ones(sz,1); %When we detect an artifact, remove 300 datapoints prior (build up of artifact) and 2000 datapoints after (3.5 sec)
        index = index+2000;
    else
        index = index +1;
    end
end

%Filter out artifacts & replace with the mean of the channels' medians
%hpc_mean = mean(median([pyra, abov, belo]));
%pfc_mean = mean(median([shal, deep]));

%pyra(outliers,:) = hpc_mean;
%abov(outliers,:) = hpc_mean;
%belo(outliers,:) = hpc_mean;
%shal(outliers,:) = pfc_mean;
%deep(outliers,:) = pfc_mean;

%% Save stuff

%result_file = fullfile(results_folder, "HPC_"+animal+"_bel_clean.continuous.mat");
%save(result_file, 'belo')
%result_file = fullfile(results_folder, "HPC_"+animal+"_pyr_clean.continuous.mat");
%save(result_file, 'pyra')