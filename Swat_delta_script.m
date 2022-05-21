%% I. Load the data from all the rats, select the top channel
%Load the libraries
addpath(genpath('C:\Users\students\Documents\Swatantra\FMAToolbox-master'));  % FMAT Toolbox

%Data path
data_path = "C:\Users\students\Documents\Swatantra\PFC_all_animals";
data_path_hpc="C:\Users\students\Documents\Swatantra\HPC_all_animals";
%Data saving location
save_path = 'C:\Users\students\Documents\Swatantra\Delta_Detections';

%Select the top channel
ratIDs = [2 3 4 5 9 10 11 201 203 204 205 206 207 209 210 211 212 213 214];  %Rat IDs
channel_ids = [ 37 37 37 37 8 8 8 7 7 7 7 7 7 7 7 7 7 7 7];%superficial channel ids
channel_ids_deep=[60 60 60 60 20 20 20 19 19 19 19 19 19 19 19 19 19 19 19]; % channels in the deep layer
ALIGNtot = [26 19 26 10 277 206 247 54 10 11 10 3 16 9 6 6 5 1 0]; % Time in minutes with respect to rat 214.
Veh = [0 1 1 0 1 0 0 1 1 0 0 1 0 0 1 1 0 1 0];  % Vector logical IDs
Cbd = [Veh == 0];   % Treatment logical IDs

artifact_thr=[3450 4000 5000 4500 4500 4500 4500 5000 5500 5000 5500 5000 6000 4500 4500 4500 4500 4500 4500 ];  % Artifact thresholds per rat provided by Victor

%Load data
% note: HPC data is loaded to look for artifact thresholds

PFC = {};
cd(data_path);  % move to the path where the data is
for i = 1:length(ratIDs)
    clear PFC_tmp
    clear HPC_tmp
    
    %Superficial PFC
    PFC_path = data_path + "\PFC_" + ratIDs(i) + "_CH"+ channel_ids(i) +".continuous.mat";    %channel 5 from end
    PFC_tmp{1} = importdata(PFC_path);

    %Deep PFC
    PFC_path = data_path + "\PFC_" + ratIDs(i) + "_CH"+ channel_ids_deep(i) +".continuous.mat";    %channel 5 from end
    PFC_tmp{2} = importdata(PFC_path);
    cd(data_path_hpc);
    
    %Above
    HPC_path = data_path_hpc + "\HPC_" + ratIDs(i) + "_CH"+ "ab.continuous.mat";    %channel 5 from end
    HPC_tmp{1} = importdata(HPC_path);

    %Pyr
    HPC_path = data_path_hpc + "\HPC_" + ratIDs(i) + "_CH"+ "pyr.continuous.mat";    %channel 5 from end
    HPC_tmp{2} = importdata(HPC_path);
    
    %Below
    HPC_path = data_path_hpc + "\HPC_" + ratIDs(i) + "_CH"+ "bel.continuous.mat";    %channel 5 from end
    HPC_tmp{3} = importdata(HPC_path);
    
  
   %Correct channels' length mismatch
    sz = min(cellfun(@length, [HPC_tmp PFC_tmp]));
    HPC_tmp{1} =HPC_tmp{1}(1:sz);
    HPC_tmp{2} =HPC_tmp{2}(1:sz);
    HPC_tmp{3} =HPC_tmp{3}(1:sz);
    PFC_tmp{1} =PFC_tmp{1}(1:sz);
    PFC_tmp{2} =PFC_tmp{2}(1:sz);
    
    
    TOT=abs(HPC_tmp{1})+abs(HPC_tmp{2})+abs(HPC_tmp{3})+abs(PFC_tmp{1})+abs(PFC_tmp{2});
    L = length(TOT);
    
    %Artifact removal
    tr = artifact_thr(i); %Visual threshold 
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
   PFC_tmp=cell2mat(PFC_tmp);
   PFC_tmp=PFC_tmp(:,1);
   PFC_tmp(outliers,:) = mean(PFC_tmp);    
    
  
  PFC{i} = PFC_tmp((600*60*15):end); % discard first 15 minutes of the data
    cd(data_path);
end
%% II.  Loop through each signal, detect delta and save in a temp variable
%initiate the return values
delta = {};
signal_filt = {};



for i = 1:length(ratIDs)
    
    % filter the signal
    fn=600; % signal frequency
    Wn1=[1/(fn/2) 6/(fn/2)]; % Cutoff=320 Hz
    [a,b] = butter(3,Wn1); %Filter coefficients
    signal_tmp = filtfilt(a,b,PFC{i}); %filtered signal
    signal_filt{i} = [zeros(ALIGNtot(i)*60*600,1); signal_tmp];  %allignment correction and saving the signal
    
    % Create the timestamp vector
    timevector = (1:length(signal_tmp))/600;% timestamps
    timevector = timevector'; %seconds
    
    % signal for detection, with the timestamps
    signal = [timevector,signal_tmp];
    
     
    % Detect delta
    delta_tmp = FindDeltaWaves(signal,i);
    
    %Alignment correction
    delta_tmp(:,1:3)=delta_tmp(:,1:3)+ALIGNtot(i)*60; 
    
    %save delta
    delta{i} = delta_tmp;
    
end
%% III. Save all the detection in a table
%cd(save_path)
%save("Delta_all_animals_final.mat",'delta')
%% IV. Plot the events to check the quality of detection.
close all
plot(signal(:,1),signal(:,2))
hold on
stem(delta_tmp(:,2),ones(size(delta_tmp(:,2)))*1000)
stem(delta_tmp(:,1),ones(size(delta_tmp(:,1)))*1000)
stem(delta_tmp(:,3),ones(size(delta_tmp(:,3)))*1000)

% %%
% %ratIDs = [2 3 4 5 9 10 11 201 203 204 205 206 207 209 210 211 212 213 214]
% %rat    = [1 2 3 4 5  6  7  8   9   10  11  12  13  14  15  16  17  18  19]
close all
rat=16;

plot((1:length(signal_filt{rat}))/600,signal_filt{rat});
hold on

stem(delta{rat}(:,2),ones(size(delta{rat}(:,2)))*1000)
stem(delta{rat}(:,1),ones(size(delta{rat}(:,1)))*1000)
stem(delta{rat}(:,3),ones(size(delta{rat}(:,3)))*1000)
% 
% % xlim([delta{rat}(16,1) delta{rat}(16,3)])


%% Binning
C=[];
vec_counts = [];
cbd_counts = [];


for i=1:length(ratIDs)
    
    C(i,:)=histcounts(delta{i}(:,2)/60,[0:45:45*12]);
    
end   