%% ||||| SWRS TYPES |||||
% Script designed to run the analysis of Ripples, Sharp Waves and Sharp
% wave ripples
clear
clc

% SET THE ANIMAL NUMBER HERE2220
% __________
animal = "214";
batch = "2";
artifact_threshold = 4500;
ripple_threshold = 5;
sharpwave_threshold = 3.5;
% ----------

source_path = 'D:\Dev\MATLAB\GenzelLab'; %uigetdir('','Source Path');

% Sleep Oscilations Tools (self)
addpath(fullfile(source_path, "sleep_oscil_tools"));

%CorticoHippocampal
addpath(genpath(fullfile(source_path, "CorticoHippocampal")))

results_folder = fullfile(source_path, "sleep_oscil_tools", "results", "detections");
if ~exist(results_folder, 'dir')
    mkdir(results_folder)
end

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

% Artifacts
artifact_thres = artifact_threshold; %Visual threshold 
outliers = false(L,1);
index = 1;
while index<L
    if TOT(index) >= artifact_thres
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
hpc_mean = mean(median([pyra, abov, belo]));
pfc_mean = mean(median([shal, deep]));

%% Sharp waves detection
% For the detection of SW (sharp waves) we use the BPL (Below Pyramidal Layer) of the HPC

%Band pass filter to better see the sharp waves
[c,d] = butter(3, [2/300 20/300]);
filtbpl = filtfilt(c, d, belo);
%Remove artifacts for detection
filtbpl(outliers) = hpc_mean;
% Detect sharp waves as signal lower than ?-5*SD
sw_thres = sharpwave_threshold;
spw = double(filtbpl <= mean(filtbpl)-sw_thres*std(filtbpl));
% Use absolute derivative for practical purposes: every sharp wave will hence start and end with a 1
dspw = abs(diff(spw));
% Get array of indices for starts and ends of sharp-waves
sw_event_marks = find(dspw);
% Get pairs of indices for starts and ends of sharp-waves
sw_event_lims = [sw_event_marks(1:2:end), sw_event_marks(2:2:end)];
% Get sharp-waves peaks
sw_num = length(sw_event_marks)/2;
peaks = zeros(sw_num,1);
for i = 1 : sw_num
    sw_lim = sw_event_lims(i,:);
    [~,I] = min(belo(sw_lim(1):sw_lim(2)));
    peaks(i) = sw_lim(1) + I -1;
end
% Get sharp-wave indexes for start, peak, end (spe)
sw_spe = [sw_event_lims(:,1) peaks sw_event_lims(:,2)];

% Detect sharp waves in the HPC layer below pyramidal layer (BPL)

%Get the local peak of each sharp wave 
%spwpeak = nan(L,1);
%index = 1;
%while index < L
%    if dspw(index) == 1
%        index2 = index;
%        index = index+1;
%        while dspw(index) == 0
%            index = index + 1;
%        end
%        locmin = min(HPC203((index2+1):(index+1),3));
%        indmin = find(HPC203(:,3) == locmin);
%        spwpeak(indmin) = locmin;
%        index = index+1;
%    else
%        index = index+1;
%    end
%end

% Sharp Waves characterization
%N_SW = sum(spwpeak); %Number of SW
%sw_epoch = ~isnan(spwpeak);
%sw_epoch = reshape(sw_epoch(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))); %separate in 1min bins
%sw_fhz = sum(sw_epoch); %Nb of SW per second for each 1min bin
%figure
%plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),length(sw_fhz)), sw_fhz)
%title('Sharp Wave rate (/sec) per 1 min bin')
% What we observe: bursts in SW activity for every NREM episode, with an
% apparent decrease as the bout lasts

%% Ripples detection
%For ripple detection, we use the pyramidal layer of the HPC
%First we bandpass for ripple frequency range
[e,f] = butter(3, [90/300 200/300]);
ripplespyr = filtfilt(e,f,pyra);
%Remove artifacts for detection
ripplespyr(outliers) = hpc_mean;
yourtimevector = (1:length(ripplespyr))/600;
r_thres = ripple_threshold;
thresh = mean(ripplespyr) + r_thres*std(ripplespyr); %threshold for ripple detection
[S, E, M] = findRipplesLisa(ripplespyr', yourtimevector, thresh, (thresh)*(1/2), 600 ); %Adrian's detection
r_spe = [S'*fn, M'*fn, E'*fn];
r_num = length(M);

%% Oscil Table creation
sw_mark = ones(sw_num,1);
r_mark = ones(r_num,1)+1;
% Columns of the table
Type = cat(1,sw_mark,r_mark);
% Start with respect to the 5 second window
Start = cat(1, sw_spe(:,1), r_spe(:,1));
%Peak with respect to raw signal
Peak = cat(1, sw_spe(:,2), r_spe(:,2));
% End with respect to the 5 second window
End = cat(1, sw_spe(:,3), r_spe(:,3));
% Five Second window with respect to raw signal
Six_Start = Peak - 1800;
Six_Start = Start - Six_Start;
Six_End = Peak + 1800;
Six_End = 3601 - (Six_End - End);
% Sleep States
Sleep_State = zeros(length(Peak),1);
for i = 1:length(Peak) 
    j = int32(Peak(i));
    Sleep_State(i) = sleep_states(j);
end
% Table
oscil_table = table(Type, Start, Peak, End, Six_Start, Six_End, Sleep_State);
% Sort by the peak of the events
oscil_table = sortrows(oscil_table, 3);

% Find last event (las row) where detection belongs to the first 15 minutes
% of the signal
start_limit = find(oscil_table{:,3}<=15*60*fn, 1, 'last');
oscil_table((1:start_limit),:) = [];

%% Find overlapping events

% Fix event arrays indices
r = oscil_table(:, [1 2 4]);
sw = oscil_table(:, [1 2 4]);
for i = 1:height(oscil_table)
    if oscil_table.Type(i) == 1
        r(i, 2) = {-2};
        r(i, 3) = {-1};
    elseif oscil_table.Type(i) == 2
        sw(i, 2) = {-4};
        sw(i, 3) = {-3};
    end
end
% Ripples
ripple = r{:,[2 3]};
% Sharp-waves
sharpwave = sw{:,[2 3]};

complexes = find_overlap(ripple, sharpwave);

for i = 1:length(complexes)
    complex = complexes{i};
    if length(complex) == 1
        for j = 1:length(complex)
            oscil_table(complex(j), 1) = {3};
        end
        oscil_table(i, 1) = {4};
    elseif length(complex) > 1
        for j = 1:length(complex)
            oscil_table(complex(j), 1) = {5};
        end 
        oscil_table(i, 1) = {6};
    end 
end

%% Waveforms 

space = zeros(fn*6+1, height(oscil_table));
wave_forms = struct('HPCabov', space, 'HPCpyra', space, 'HPCbelo', space, 'PFCshal', space, 'PFCdeep', space);
for i = 1: height(oscil_table)
    peak = oscil_table{i,3};
    if peak-1800 < 1
        lolim = 1;
        hilim = 3600;
    elseif peak+1800 > L
        lolim = L-3600; 
        hilim = L;
    else
        lolim = peak-1800; 
        hilim = peak+1800;
    end
    window = fix([lolim, hilim]);
    wave_forms.HPCpyra(:,i) = pyra(window(1) : window(2));
    wave_forms.HPCabov(:,i) = abov(window(1) : window(2));
    wave_forms.HPCbelo(:,i) = belo(window(1) : window(2));
    wave_forms.PFCshal(:,i) = shal(window(1) : window(2));
    wave_forms.PFCdeep(:,i) = deep(window(1) : window(2));
end

%% Grouped table

% Singles pretable
singles = oscil_table(oscil_table.Type == 1 | ...
                   oscil_table.Type == 2,:);
singles.Properties.VariableNames(1) = {'Form'};
singles.Form = num2cell(singles.Form);
Type = repelem("singles",height(singles));
singles = [table(Type') singles];
singles.Properties.VariableNames(1) = {'Type'};

% Simples and complex pretables
Type = string(zeros(1,length(complexes)));
Form = cell(1,length(complexes));
Start = zeros(1,length(complexes));
Peak = zeros(1,length(complexes));
End = zeros(1,length(complexes));
Six_Start = zeros(1,length(complexes));
Six_End = zeros(1,length(complexes));
for i = 1:length(complexes)
    complex = complexes{i};
    type = length(complex);
    if type == 1 || type > 1
        indexes = sort([i complex]);
        if type > 1; Type(i) = "complex"; else; Type(i) = "simple"; end
        Form{i} = oscil_table{indexes,1}';
        Start(i) = min(oscil_table{indexes,2});
        Peak(i) = oscil_table{i,3};
        End(i) = max(oscil_table{indexes,4});
        Six_Start(i) = oscil_table{i,5};
        Six_End(i) = oscil_table{i,6};
    end 
end
indexes = find(~strcmp("0", Type));
Type = Type(indexes)';
Form = Form(indexes)';
Start = Start(indexes)';
Peak = Peak(indexes)';
End = End(indexes)';
Six_Start = Six_Start(indexes)';
Six_End = Six_End(indexes)';

Sleep_State = zeros(length(Peak),1);
for i = 1:length(Peak) 
    j = int32(Peak(i));
    Sleep_State(i) = sleep_states(j);
end

simples_complexes = table(Type, Form, Start, Peak, End, Six_Start, Six_End, Sleep_State);

% Complex table
grouped_oscil_table = [singles; simples_complexes];
% Sort by the peak of the events
grouped_oscil_table = sortrows(grouped_oscil_table, 4);

%% Grouped Waveforms 

space = zeros(fn*6+1, height(grouped_oscil_table));
grouped_wave_forms = struct('HPCabov', space, 'HPCpyra', space, 'HPCbelo', space, 'PFCshal', space, 'PFCdeep', space);
for i = 1: height(grouped_oscil_table)
    peak = grouped_oscil_table{i,4};
    if peak-1800 < 1
        lolim = 1;
        hilim = 3600;
    elseif peak+1800 > L
        lolim = L-3600; 
        hilim = L;
    else
        lolim = peak-1800; 
        hilim = peak+1800;
    end
    window = fix([lolim, hilim]);
    grouped_wave_forms.HPCpyra(:,i) = pyra(window(1) : window(2));
    grouped_wave_forms.HPCabov(:,i) = abov(window(1) : window(2));
    grouped_wave_forms.HPCbelo(:,i) = belo(window(1) : window(2));
    grouped_wave_forms.PFCshal(:,i) = shal(window(1) : window(2));
    grouped_wave_forms.PFCdeep(:,i) = deep(window(1) : window(2));
end

%% Visualization of events

% Sharp Waves
sw_evs = oscil_table(oscil_table.Type == 1,:);
sw_evs = (sw_evs{:,[2 4]}/600)';
sw_pos = zeros(size(sw_evs));
sw_peaks = oscil_table(oscil_table.Type == 1,:);
sw_peaks = sw_peaks{:,3}/600;
M_dur_sw_p = NaT(1, length(sw_peaks)) - NaT(1);
M_dur_sw = NaT(2, length(sw_evs)) - NaT(1);
for kk=1:length(sw_evs)
    M_dur_sw_p(kk) = duration([0 0 sw_peaks(kk)]);
    M_dur_sw(1,kk) = duration([0 0 sw_evs(1,kk)]);
    M_dur_sw(2,kk) = duration([0 0 sw_evs(2,kk)]);
end

% Ripples
r_evs = oscil_table(oscil_table.Type == 2,:);
r_evs = (r_evs{:,[2 4]}/600)';
r_pos = zeros(size(r_evs))+1;
r_peaks = oscil_table(oscil_table.Type == 2,:);
r_peaks = r_peaks{:,3}/600;
M_dur_r_p = NaT(1, length(r_peaks)) - NaT(1);
M_dur_r = NaT(2, length(r_evs)) - NaT(1);
for kk=1:length(r_evs)
    M_dur_r_p(kk) = duration([0 0 r_peaks(kk)]);
    M_dur_r(1,kk) = duration([0 0 r_evs(1,kk)]);
    M_dur_r(2,kk) = duration([0 0 r_evs(2,kk)]);
end

% SW_SWR
sw_swr_evs = oscil_table(oscil_table.Type == 3 | oscil_table.Type == 5,:);
sw_swr_evs = (sw_swr_evs{:,[2 4]}/600)';
sw_swr_pos = zeros(size(sw_swr_evs))+3;
sw_swr_peaks = oscil_table(oscil_table.Type == 3 | oscil_table.Type == 5,:);
sw_swr_peaks = sw_swr_peaks{:,3}/600;
M_dur_sw_swr_p = NaT(1, length(sw_swr_peaks)) - NaT(1);
M_dur_sw_swr = NaT(2, length(sw_swr_evs)) - NaT(1);
for kk=1:length(sw_swr_evs)
    M_dur_sw_swr_p(kk) = duration([0 0 sw_swr_peaks(kk)]);
    M_dur_sw_swr(1, kk) = duration([0 0 sw_swr_evs(1,kk)]);
    M_dur_sw_swr(2, kk) = duration([0 0 sw_swr_evs(2,kk)]);
end

% R_SWR
r_swr_evs = oscil_table(oscil_table.Type == 4 | oscil_table.Type == 6,:);
r_swr_evs = (r_swr_evs{:,[2 4]}/600)';
r_swr_pos = zeros(size(r_swr_evs))+2;
r_swr_peaks = oscil_table(oscil_table.Type == 4 | oscil_table.Type == 6,:);
r_swr_peaks = r_swr_peaks{:,3}/600;
M_dur_r_swr_p = NaT(1, length(r_swr_peaks)) - NaT(1);
M_dur_r_swr = NaT(2, length(r_swr_evs)) - NaT(1);
for kk=1:length(r_swr_evs)
    M_dur_r_swr_p(kk) = duration([0 0 r_swr_peaks(kk)]);
    M_dur_r_swr(1,kk) = duration([0 0 r_swr_evs(1,kk)]);
    M_dur_r_swr(2,kk) = duration([0 0 r_swr_evs(2,kk)]);
end

% Display all prefiltered signals
figure
tiledlayout(5,1)
increment = 3;

tt1 = nexttile;
sig = increment.*abov;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), sig, 'color', [0.3010, 0.7450, 0.9330])
title('HPC - above pyramidal layer')
minsig = min(sig);
maxsig = max(sig);
ylim([minsig maxsig+1200])
hold on
plot(M_dur_r, (maxsig)*ones(size(M_dur_r)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw, (maxsig+750)*ones(size(M_dur_sw)), 'color', [1, 0, 0], 'LineWidth', 12)
hold on
plot(M_dur_r_swr, (maxsig)*ones(size(M_dur_r_swr)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw_swr, (maxsig+750)*ones(size(M_dur_sw_swr)), 'color', [1, 0, 0], 'LineWidth', 12)

tt2 = nexttile;
sig = increment.*pyra;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), sig, 'color', [0.3010, 0.7450, 0.9330])
title('HPC - pyramidal layer')
hold on
r_p = stem(M_dur_r_p,500*ones(size(M_dur_r_p)), 'color', [0.4660, 0.6740, 0.1880]);
set(r_p, 'Marker', 'none')
hold on
r_swr_p = stem(M_dur_r_swr_p,500*ones(size(M_dur_r_swr_p)), 'color', [0.4660, 0.6740, 0.1880]);
set(r_swr_p, 'Marker', 'none')
hold on
minsig = min(sig);
maxsig = max(sig);
ylim([minsig maxsig+1200])
hold on
plot(M_dur_r, (maxsig)*ones(size(M_dur_r)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw, (maxsig+750)*ones(size(M_dur_sw)), 'color', [1, 0, 0], 'LineWidth', 12)
hold on
plot(M_dur_r_swr, (maxsig)*ones(size(M_dur_r_swr)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw_swr, (maxsig+750)*ones(size(M_dur_sw_swr)), 'color', [1, 0, 0], 'LineWidth', 12)


tt3 = nexttile;
sig = increment.*belo;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), sig, 'color', [0.3010, 0.7450, 0.9330])
title('HPC - below pyramidal layer')
hold on
sw_p = stem(M_dur_sw_p,300*ones(size(M_dur_sw_p)), 'color', [1, 0, 0]);
set(sw_p, 'Marker', 'none')
hold on
sw_swr_p = stem(M_dur_sw_swr_p,300*ones(size(M_dur_sw_swr_p)), 'color', [1, 0, 0]);
set(sw_swr_p, 'Marker', 'none')
hold on
minsig = min(sig);
maxsig = max(sig);
ylim([minsig maxsig+3000])
hold on
plot(M_dur_r, (maxsig+200)*ones(size(M_dur_r)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw, (maxsig+1850)*ones(size(M_dur_sw)), 'color', [1, 0, 0], 'LineWidth', 12)
hold on
plot(M_dur_r_swr, (maxsig+200)*ones(size(M_dur_r_swr)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw_swr, (maxsig+1850)*ones(size(M_dur_sw_swr)), 'color', [1, 0, 0], 'LineWidth', 12)

tt4 = nexttile;
sig = increment.*shal;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), sig, 'color', [0.3010, 0.7450, 0.9330])
title('PFC - shallow layer')
hold on
minsig = min(sig);
maxsig = max(sig);
ylim([minsig maxsig+1200])
hold on
plot(M_dur_r, (maxsig)*ones(size(M_dur_r)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw, (maxsig+750)*ones(size(M_dur_sw)), 'color', [1, 0, 0], 'LineWidth', 12)
hold on
plot(M_dur_r_swr, (maxsig)*ones(size(M_dur_r_swr)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw_swr, (maxsig+750)*ones(size(M_dur_sw_swr)), 'color', [1, 0, 0], 'LineWidth', 12)


tt5 = nexttile;
sig = increment.*deep;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), sig, 'color', [0.3010, 0.7450, 0.9330])
title('PFC - deep layer')
hold on
hold on
plot(M_dur_r, (maxsig)*ones(size(M_dur_r)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw, (maxsig+750)*ones(size(M_dur_sw)), 'color', [1, 0, 0], 'LineWidth', 12)
hold on
plot(M_dur_r_swr, (maxsig)*ones(size(M_dur_r_swr)), 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 12)
hold on
plot(M_dur_sw_swr, (maxsig+750)*ones(size(M_dur_sw_swr)), 'color', [1, 0, 0], 'LineWidth', 12)

linkaxes([tt1 tt2 tt3 tt4 tt5], 'x')

%visualization_path = fullfile(results_path, "visualization" + animal + ".fig");
%savefig(visualization_path);

%% Save everything
dataset_filename = fullfile(results_folder,"detections" + animal + ".mat");
save(dataset_filename, 'oscil_table', 'wave_forms', 'grouped_oscil_table', 'grouped_wave_forms')

%%
% We have the start/end time of each ripple, we create one long logical
% vector where 1's are our ripples
% ripple_hpc203 = zeros(1,L);
% 
% index = 1; % Timestamp index.
% index2 = 1; %Events index.
% while index < length(yourtimevector) & index2<length(E)
%     if S(index2) == yourtimevector(index)
%         while E(index2) ~= yourtimevector(index)
%             ripple_hpc203(index) = 1;
%             index = index + 1;
%         end
%         index2 = index2 + 1;
%     else
%         index = index + 1;
%     end
% end
% 
% % Now we create a "Venn diagram" of ripples alone, SW alone and co-occuring
% % SWripples (or SWR)
% % We group by bins of 200ms
% wl = 600*0.2; %200ms bins ==> 600 points = 1sec
% tmp = floor(L/wl);
% 
% bin_PL = reshape(HPC203(1:tmp*wl,1), wl, tmp); %Pyramidal Layer
% bin_BPL = reshape(HPC203(1:tmp*wl,3), wl, tmp); % Below Pyramidal Layer
% bin_ripples = reshape(ripple_hpc203(1:tmp*wl), wl, tmp); 
% bin_filt_rip = reshape(ripplespyr(1:tmp*wl), wl, tmp); % Filtered pyramidal layer for ripples
% tmpsw = ~isnan(spwpeak);
% bin_sw = reshape(tmpsw(1:tmp*wl), wl, tmp);
% 
% % Sum per bin ==> if sum is >0, means that one event is happening (either
% % SW in the BPL or ripple in the pyramidal layer - PL)
% sumbin_sw = sum(bin_sw); 
% sumbin_rip = sum(bin_ripples);
% 
% venn = zeros(1,tmp);
% for kk=1:tmp
%     venn(kk) = 1*(sumbin_sw(kk)>0 & sumbin_rip(kk) == 0) + 2*(sumbin_sw(kk) == 0 & sumbin_rip(kk) > 0) + 3*(sumbin_sw(kk) > 0 & sumbin_rip(kk) > 0);
% end
% % So in the end that is: 1 for SWs alone, 2 for ripples alone, 3 for SWR
% 
% 
% COUNT = zeros(1,3);
% 
% % Need to adjust venn for consecutive "events" = one longer event
% index1 = 1;
% while index1 < length(venn)
%     if venn(index1) ~=0 % Beginning of an "event"
%        index2 = index1+1;
%        while venn(index2) ~=0 & index2 < length(venn) %We look the length of that event (or consecutive events)
%            index2 = index2+1;
%        end
%        event = venn(index1:(index2 -1));
%        if length(event) >= 2 %If that event last longer than 2 200ms-bins, then it potentially is a mix of different "types" of events
%            % Several cases
%            if ismember(3, event)
%                venn(index1:(index2-1)) = 3; %If one 200 ms bin is a SWR, then the whole event is a SWR
%            elseif ismember(1, event) & ismember(2, event) 
%                 venn(index1:(index2-1)) = 3; %If it contains both a RIP without SW and a SW without RIP, then it's a SWR
%            end
%        end
%        COUNT = COUNT + ismember(1,venn(index1:(index2-1)))*[1 0 0] + ismember(2,venn(index1:(index2-1)))*[0 1 0]+ ismember(3,venn(index1:(index2-1)))*[0 0 1];
%        index1 = index2;
%     else
%         index1 = index1 +1;
%     end
% end
% 
% % Barplot of the number of events
% figure
% b = bar([0.5 2 3.5], COUNT);
% b.FaceColor = 'flat';
% b.EdgeColor = 'none';
% xticks([0.5 2 3.5])
% xticklabels({'SW without ripple','Ripple without SW', 'SWR'})
% xtickangle(25)
% b.CData(1,:) = [0    0.4471    0.7412];
% b.CData(2,:) = [0.9216    0.2863    0.2863];
% b.CData(3,:) = [0.4667    0.7804    0.3529];
% title('Number of each event type for rat 203')
% 
% % Plots
% % Note - Patches can be really long to display in plots
% figure
% tiledlayout(4,1)
% % t1 = nexttile;
% colorsvenn = ['b', 'r', 'g'];
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), HPC203(:,1))
% % hold on
% % for ii=1:length(venn)
% %     if venn(ii) ~=0
% %         patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(HPC203(:,1))) -max(abs(HPC203(:,1))) max(abs(HPC203(:,1))) max(abs(HPC203(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
% %         hold on
% %     end
% % end
% % title('Pyramidal Layer')
% t1 = nexttile;
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), ripplespyr)
% hold on
% stem(M_dur,600*ones(size(M_dur)))
% title('Pyramidal layer filtered for ripples')
% t2 = nexttile;
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), HPC203(:,3))
% hold on
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), filtblwpyr)
% % hold on
% scatter(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L),spwpeak, 'red', 'filled')
% title('Below pyramidal layer filtered with sharp waves')
% t3 = nexttile;
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC203(:,1))
% hold on
% for ii=1:length(venn)
%     if venn(ii) ~=0
%         patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(HPC203(:,1))) -max(abs(HPC203(:,1))) max(abs(HPC203(:,1))) max(abs(HPC203(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
%         hold on
%     end
% end
% title('PFC channel')
% t4 = nexttile;
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), [0; diff(PFC203(:,1))])
% linkaxes([t1 t2 t3 t4], 'x')
% 
% 
% %% Amplitudes of sharp waves vs ripples
% % filtbpl for the sharp waves
% % ripplespyr for the ripples
% amplitudes = [];
% for ii=1:length(venn)
%    if venn(ii)~=0 
%        sw_vec = filtbpl(((ii-1)*120+1):ii*120); % 1 datapoint in venn = 120 datapoints in the raw signal (200ms * 600 datapoints/sec)
%        rip_vec = ripplespyr(((ii-1)*120+1):ii*120);
%        pfc1 = PFC203(((ii-1)*120+1):ii*120,1);
%        pfc2 = PFC203(((ii-1)*120+1):ii*120,2);
%        ampl_sw = mean(filtbpl)-min(sw_vec);
%        ampl_rip = peak2peak(rip_vec);
%        ampl_pfc1 = peak2peak(pfc1);
%        ampl_pfc2 = peak2peak(pfc2);
%        amplitudes = [amplitudes; [ampl_sw ampl_rip ampl_pfc1 ampl_pfc2]];
%    end
% end
% 
% figure
% scatter(amplitudes(:,1), amplitudes(:,2))
% title('Amplitude of sharp waves vs ripples')
% xlabel('Sharp wave amplitude') 
% ylabel('Ripple amplitude')
% 
% %No clear tendency, although there seems to be some kind of correlation
% %% REMINDER
% % Venn: 1=SW, blue, 2=rip, red, 3=SWR, green
% %% Spectral power in PFC after events 
% 
% % Spectrogram, using an overlap window of 1sec (timestep 0.2 sec)
% params=struct('tapers',[],'Fs',[],'fpass',[]);
% params.tapers = [5/2 4];
% params.Fs = 600;
% params.fpass = [1 200];
% [S,t,f] = mtspecgramc(PFC203(:,1),[1 0.2],params);
% 
% figure
% tiledlayout(2,1)
% t1 = nexttile;
% surf(duration(0,0,t),f,10*log10(S'),'EdgeColor', 'none')
% view(2)
% colormap(t1, 'parula')
% %colorbar();
% set(gca,'YDir','normal')
% title('Spectrogram')
% 
% t2 = nexttile;
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC203(:,1))
% hold on
% % Color patch very long to plot, find a better way to plot it
% % for ii=1:length(venn)
% %     if venn(ii) ~=0
% %         patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(PFC203(:,1))) -max(abs(PFC203(:,1))) max(abs(PFC203(:,1))) max(abs(PFC203(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
% %         hold on
% %     end
% % end
% title('PFC channel superficial')
% 
% linkaxes([t1 t2], 'x')
% 
% 
% %% Average over a 1sec (?) window for every event type
% % We align on each event, and plot the average PFC spectrogram 2 sec around the
% % event
% window_size = 2; %sec
% 
% av_pfc_power_rip = zeros(window_size*5,340);
% av_pfc_power_sw = zeros(window_size*5,340);
% av_pfc_power_swr = zeros(window_size*5,340);
% n_sw = 0;
% n_rip = 0;
% n_swr = 0;
% 
% index1 = 1;
% while index1 < length(venn)
%     if venn(index1) ~=0 % Beginning of an "event"
%        index2 = index1+1;
%        while venn(index2) ~=0 & index2 < length(venn) %We look at the length of that event (or consecutive events)
%            index2 = index2+1;
%        end
%        event = venn(index1:(index2 -1));
%        dur = (index2-index1)*0.2;
%        delta = window_size-dur;
%        [~,window_low] = min(abs(t-(index1*0.2-delta/2))); %Find index in t corresponding to the timing of that event 
%        [~,window_high] = min(abs(t-((index2-1)*0.2 + delta/2)));
%        pfc_spectrogram_window = S(window_low:window_high,:);
%        if event(1) == 1;
%            av_pfc_power_sw = av_pfc_power_sw + pfc_spectrogram_window;
%            n_sw = n_sw + 1;
%        elseif event(1) == 2;
%            av_pfc_power_rip = av_pfc_power_rip + pfc_spectrogram_window;
%            n_rip = n_rip + 1;
%        elseif event(1) == 3;
%            av_pfc_power_swr = av_pfc_power_swr + pfc_spectrogram_window;
%            n_swr = n_swr + 1;
%        end
%        index1 = index2;
%     else
%         index1 = index1 +1;
%     end
% end
% av_pfc_power_sw = av_pfc_power_sw / n_sw;
% av_pfc_power_rip = av_pfc_power_rip / n_rip;
% av_pfc_power_swr = av_pfc_power_swr / n_swr;
% 
% figure
% tiledlayout(1,3)
% t1 = nexttile();
%     surf(duration(0,0,(1:window_size*5)*0.2),f, log(av_pfc_power_rip'), 'EdgeColor','none')
%     view(2)
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around ripples')
% t2 = nexttile();
%      surf(duration(0,0,(1:window_size*5)*0.2),f, log(av_pfc_power_sw'), 'EdgeColor','none')
%     view(2)
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around SWs')
% t3 = nexttile();
%      surf(duration(0,0,(1:window_size*5)*0.2),f, log(av_pfc_power_swr'), 'EdgeColor','none')
%     view(2)
%     colormap(t1, 'parula')
%     %colorbar();
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around SWRs')
%     
%     linkaxes([t1 t2 t3], 'y')
%     
% %% Average over a 1sec window for every event type - split low/high fhz
% % Same plot but separate low/high frequencies
% %1-20 Hz
% 
% [~,fmin] = min(abs(f-1));
% [~,fmax] = min(abs(f-20));
% figure
% tiledlayout(1,3)
% t1 = nexttile();
%     surf(duration(0,0,(1:window_size*5)*0.2),f(fmin:fmax), log(av_pfc_power_rip(:,fmin:fmax)'), 'EdgeColor','none')
%     view(2)
%     ylim([0 20])
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around ripples')
% t2 = nexttile();
%      surf(duration(0,0,(1:window_size*5)*0.2),f(fmin:fmax), log(av_pfc_power_sw(:,fmin:fmax)'), 'EdgeColor','none')
%     view(2)
%     ylim([0 20])
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around SWs')
% t3 = nexttile();
%      surf(duration(0,0,(1:window_size*5)*0.2),f(fmin:fmax), log(av_pfc_power_swr(:,fmin:fmax)'), 'EdgeColor','none')
%     view(2)
%     ylim([0 20])
%     colormap(t1, 'parula')
%     %colorbar();
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around SWRs (0-20Hz)')
%     
%     linkaxes([t1 t2 t3], 'y')
%     
%   
% % 20 - 200 Hz
% 
% figure
% tiledlayout(1,3)
% t1 = nexttile();
%     surf(duration(0,0,(1:window_size*5)*0.2),f(fmax:end), log(av_pfc_power_rip(:,fmax:end)'), 'EdgeColor','none')
%     view(2)
%     ylim([20 200])
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around ripples')
% t2 = nexttile();
%      surf(duration(0,0,(1:window_size*5)*0.2),f(fmax:end), log(av_pfc_power_sw(:,fmax:end)'), 'EdgeColor','none')
%     view(2)
%     ylim([20 200])
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around SWs')
% t3 = nexttile();
%      surf(duration(0,0,(1:window_size*5)*0.2),f(fmax:end), log(av_pfc_power_swr(:,fmax:end)'), 'EdgeColor','none')
%     view(2)
%     ylim([20 200])
%     colormap(t1, 'parula')
%     %colorbar();
%     set(gca,'YDir','normal')
%     title('Averaged spectrogram of PFC activity around SWRs (20-200Hz)')
%     
%     linkaxes([t1 t2 t3], 'y')
%     
%     
% %% Spectral power in PFC after events - with the other PFC channel
% 
% % Spectrogram, using an overlap window of 1sec (timestep 0.2 sec)
% params=struct('tapers',[],'Fs',[],'fpass',[]);
% params.tapers = [5/2 4];
% params.Fs = 600;
% params.fpass = [1 200];
% [S,t,f] = mtspecgramc(PFC203(:,2),[1 0.2],params);
% 
% figure
% tiledlayout(2,1)
% t1 = nexttile;
% surf(duration(0,0,t),f,10*log10(S'),'EdgeColor', 'none')
% view(2)
% colormap(t1, 'parula')
% %colorbar();
% set(gca,'YDir','normal')
% title('Spectrogram')
% 
% t2 = nexttile;
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC203(:,1))
% hold on
% % Color patch very long to plot, find a better way to plot it
% % for ii=1:length(venn)
% %     if venn(ii) ~=0
% %         patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(PFC203(:,1))) -max(abs(PFC203(:,1))) max(abs(PFC203(:,1))) max(abs(PFC203(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
% %         hold on
% %     end
% % end
% title('PFC channel superficial')
% 
% linkaxes([t1 t2], 'x')
% 
% %% Looking at features of sw/rip
% % Here we look at every event and we retrieve information such as
% % amplitude, mean freq, durations
% venn_extended = [repelem(venn,120) zeros(1,L-length(venn)*120)]; %Extend venn to align with raw signals
% 
% % We begin with ripples characteristics
% ripples_charac = [];
% index1 = 1;
% while index1 < L
%     if ripple_hpc203(index1) == 1
%         index2 = index1;
%         while ripple_hpc203(index2) == 1 & index2 < L
%             index2 = index2 + 1;
%         end
%         AMPL = peak2peak(ripplespyr(index1:(index2-1)));
%         [pxx,f] = pmtm(ripplespyr(index1:(index2-1)),4,[],600);
%         FREQ = meanfreq(pxx, 600, [70 290]);
%         DUR = (index2-index1)*1000/600;
%         TYPE = venn_extended(index1);
%         ripples_charac = [ripples_charac; AMPL DUR FREQ TYPE];
%         index1 = index2;
%     else
%       index1 = index1 + 1;  
%     end
% end
% 
% % Amplitude of ripples
% 
% figure
% h1 = histogram(ripples_charac(ripples_charac(:,4) == 3, 1),'BinWidth', 10 ,'FaceAlpha', 0.4, 'FaceColor','g');
% hold on
% h2 = histogram(ripples_charac(ripples_charac(:,4) == 2, 1), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
% title('Distribution of the ripples amplitudes')
% xlabel('Peak to trough amplitude')
% legend({'SWRs', 'Ripples only'})
% 
% % Duration of ripples
% figure
% h1 = histogram(ripples_charac(ripples_charac(:,4) == 3, 2),'BinWidth', 10 ,'FaceAlpha', 0.4, 'FaceColor','g');
% hold on
% h2 = histogram(ripples_charac(ripples_charac(:,4) == 2, 2), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
% title('Distribution of the ripples duration')
% xlabel('Duration (ms)')
% legend({'SWRs', 'Ripples only'})
% 
% % Mean frequency of ripples
% figure
% h1 = histogram(ripples_charac(ripples_charac(:,4) == 3, 3),'BinWidth', 3 ,'FaceAlpha', 0.4, 'FaceColor','g');
% hold on
% h2 = histogram(ripples_charac(ripples_charac(:,4) == 2, 3), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
% title('Distribution of the ripples mean frequency')
% xlabel('Frequency (Hz)')
% legend({'SWRs', 'Ripples only'})
% 
% % Now SW characteristics
% SW_charac = [];
% for ii = 1:L
%     if ~isnan(spwpeak(ii)) 
%         AMPL = -spwpeak(ii);
%         TYPE = venn_extended(ii);
%         SW_charac = [SW_charac; AMPL TYPE];
%     end
% end
% % Amplitude of sharp waves
% figure
% h1 = histogram(SW_charac(SW_charac(:,2) == 3, 1),'BinWidth', 50 ,'FaceAlpha', 0.4, 'FaceColor','g');
% hold on
% h2 = histogram(SW_charac(SW_charac(:,2) == 1, 1), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','b');
% title('Distribution of the SW amplitudes')
% xlabel('SW amplitude')
% legend({'SWRs', 'SWs only'})
% 
% %% Split ripples
% % Here we want to split ripples between "slow" and "fast" ripples
% M_dur = NaT(1, length(M)) - NaT(1);
% for kk=1:length(M)
%     M_dur(kk) = duration([0 0 M(kk)]);
% end
% 
% ripple_hpc203 = zeros(1,L);
% 
% % We just reconstruct a logical vector with 0 - no ripple and 1 - ripple
% index = 1;
% index2 = 1;
% while index < length(yourtimevector) & index2<length(E)
%     if S(index2) == yourtimevector(index)
%         while E(index2) ~= yourtimevector(index)
%             ripple_hpc203(index) = 1;
%             index = index + 1;
%         end
%         index2 = index2 + 1;
%     else
%         index = index + 1;
%     end
% end
% 
% % Slow ripples
% [e,f] = butter(3, [90/300 115/300]);
% rip_slow = filtfilt(e,f,HPC203(:,1));
% yourtimevector = (1:length(rip_slow))/600;
% thresh_slow = mean(rip_slow) + 5*std(rip_slow);
% [S_slow, E_slow, M_slow] = findRipplesLisa(rip_slow', yourtimevector, thresh_slow, (thresh_slow)*(1/2), [] );
% 
% % Fast ripples
% [e,f] = butter(3, [115/300 200/300]);
% rip_fast = filtfilt(e,f,HPC203(:,1));
% thresh_fast = mean(rip_fast) + 5*std(rip_fast);
% [S_fast, E_fast, M_fast] = findRipplesLisa(rip_fast', yourtimevector, thresh_fast, (thresh_fast)*(1/2), [] );
% 
% figure
% t = tiledlayout(2,1);
% t1 = nexttile;
% plot(rip_slow)
% t2 = nexttile
% plot(rip_fast)
% linkaxes([t1 t2], 'x')
% 
% % Features of slow/fast ripples
% 
% ripple_slow_hpc203 = zeros(1,L);
% 
% index = 1;
% index2 = 1;
% while index < length(yourtimevector) & index2<length(E_slow)
%     if S_slow(index2) == yourtimevector(index)
%         while E_slow(index2) ~= yourtimevector(index)
%             ripple_slow_hpc203(index) = 1;
%             index = index + 1;
%         end
%         index2 = index2 + 1;
%     else
%         index = index + 1;
%     end
% end
% 
% ripple_fast_hpc203 = zeros(1,length(HPC203(:,3)));
% 
% index = 1;
% index2 = 1;
% while index < length(yourtimevector) & index2<length(E_fast)
%     if S_fast(index2) == yourtimevector(index)
%         while E_fast(index2) ~= yourtimevector(index)
%             ripple_fast_hpc203(index) = 1;
%             index = index + 1;
%         end
%         index2 = index2 + 1;
%     else
%         index = index + 1;
%     end
% end
% 
% ripples_slow_charac = [];
% index1 = 1
% while index1 < L
%     if ripple_slow_hpc203(index1) == 1
%         index2 = index1;
%         while ripple_slow_hpc203(index2) == 1 & index2 < L
%             index2 = index2 + 1;
%         end
%         AMPL = peak2peak(rip_slow(index1:(index2-1)));
% %         [pxx,f] = pmtm(rip_slow(index1:(index2-1)),4,[],600);
% %         FREQ = meanfreq(pxx, 600, [70 290]);
%         DUR = (index2-index1)*1000/600;
%         ripples_slow_charac = [ripples_slow_charac; AMPL DUR];
%         index1 = index2;
%     else
%       index1 = index1 + 1;  
%     end
% end
% 
% 
% ripples_fast_charac = [];
% index1 = 1
% while index1 < L
%     if ripple_fast_hpc203(index1) == 1
%         index2 = index1;
%         while ripple_fast_hpc203(index2) == 1 & index2 < L
%             index2 = index2 + 1;
%         end
%         AMPL = peak2peak(rip_fast(index1:(index2-1)));
% %         [pxx,f] = pmtm(rip_fast(index1:(index2-1)),4,[],600);
% %         FREQ = meanfreq(pxx, 600, [70 290]);
%         DUR = (index2-index1)*1000/600;
%         ripples_fast_charac = [ripples_fast_charac; AMPL DUR];
%         index1 = index2;
%     else
%       index1 = index1 + 1;  
%     end
% end
% 
% power_slow_fast = [];
% index1 = 1
% while index1 < L
%     if ripple_hpc203(index1) == 1
%         index2 = index1;
%         while ripple_hpc203(index2) == 1 & index2 < L
%             index2 = index2 + 1;
%         end
%         [pxx,f] = pmtm(ripplespyr(index1:(index2-1)),4,[],600);
%         slopo = spectral_power(80, 110, 1, pxx, f);
%         fapo = spectral_power(110, 200, 1, pxx, f);
% 
%         power_slow_fast = [power_slow_fast; slopo fapo];
%         index1 = index2;
%     else
%       index1 = index1 + 1;  
%     end
% end
% 
% %plots
% 
% 
% % Amplitude of ripples
% 
% figure
% h1 = histogram(ripples_slow_charac(:, 1),'BinWidth', 10 ,'FaceAlpha', 0.4, 'FaceColor','g');
% hold on
% h2 = histogram(ripples_fast_charac(:, 1), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
% title('Distribution of the ripples amplitudes')
% xlabel('Peak to trough amplitude')
% legend({'Slow ripples', 'Fast ripples'})
% 
% % Duration of ripples
% figure
% h1 = histogram(ripples_slow_charac(:, 2),'BinWidth', 10 ,'FaceAlpha', 0.4, 'FaceColor','g');
% hold on
% h2 = histogram(ripples_fast_charac(:, 2), 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor','r');
% title('Distribution of the ripples duration')
% xlabel('Duration (ms)')
% legend({'Slow ripples', 'Fast ripples'})
% 
% % slow power vs fast power
% 
% figure
% scatter(log(power_slow_fast(:,1)), log(power_slow_fast(:,2)))
% xlabel('Log-transformed slow power (90-115Hz)')
% ylabel('Log-transformed fast power (115-200Hz)')
% title('Slow vs fast ripple power')
% 
% % Adrian's gaussmix_slow_fast model
% % all ripples
% freq1 = ripples_charac(:,3);
% [cutoff1, sd1] = gaussmix_slow_fast(freq1);
% 
% %Rip alone
% freq2 = ripples_charac(ripples_charac(:,4) == 2,3);
% [cutoff2, sd2] = gaussmix_slow_fast(freq2);
% 
% %SWRs only
% freq3 = ripples_charac(ripples_charac(:,4) == 3,3);
% [cutoff3, sd3] = gaussmix_slow_fast(freq3);
% 
% % ==> Gaussian mixture model does not fit: it works but the cutoff is
% % obviously wrong... then should we select it by hand?
% %% Look at events during UP/DOWN
% % % We know the PFC during NREM is dominated by UP states followed by DOWN
% % % states and so on
% 
% % IGNORE THE REST OF THIS SECTION
% % [c,d] = butter(3, [0.1/300 2/300]);
% % up_down = filtfilt(c,d,PFC203(:,1));
% % % We load the scored data
% % data = load('C:\Users\students\Documents\Tugdual\PCA_Scorer_batch_2\Rat203\states2');
% % states203 = data.clusters;
% % 
% % %Use 20-60Hz to detect DOWN states
% % epdata = data2ep(PFC203(:,1), 0.2, 600); 
% % [pxx,f] = pmtm(epdata,4,[],ds_fhz);
% % 
% % power20_60 = [];
% % for jj=1:size(epdata,2)
% %       val = spectral_power(25, 180, jj, pxx, f);
% %       power20_60 = [power20_60; val];
% % end
% % ampl_power20_60 = [repelem(power20_60, 120) ;zeros(L-120*length(power20_60),1)];
% % 
% % figure
% % tiledlayout(4,1)
% % t1 = nexttile;
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC203(:,1))
% % t2 = nexttile;
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), up_down)
% % t3 = nexttile;
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), ampl_power20_60)
% % 
% % 
% % %VISUAL THRESHOLD FOR DOWN STATES : POWER > 500
% % 
% % UP_DOWN_log = 2*(ampl_power20_60<500) + 1 *(ampl_power20_60>=500);
% % t4 = nexttile;
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L),UP_DOWN_log)
% % ylim([0 3])
% % linkaxes([t1 t2 t3 t4], 'x');
% 
% %% Other UPDOWN detection
% % We detect UPDOWN states with the hilbert transform, then retrieve the
% % phase
% 
% %First filter for low frequencies
% [c,d] = butter(3, [0.1/300 1/300]);
% up_down2 = filtfilt(c,d,PFC203(:,1));
% %Get only NREM sleep, else during REM the hilbert will do weird stuff, here
% %we take by hand the longest NREM bout
% portion_nrem = up_down2(10620000:12620000);
% t=1:length(portion_nrem);
% x=portion_nrem;
% %Compute hilbert transform
% h=hilbert(x);
% %Get phase and turn into degrees
% phase_radians=angle(h);
% de=rad2deg(phase_radians);
% 
% figure
% plot(t,de)
% phase_degrees=mod(de,360);
% plot(t,x); hold on; plot(t,phase_degrees)
% title('Updown detection')
% legend('Filtered signal','Phase of hilbert transform')
% 
% % Cut the venn to fit the selectedNREM period 
% venn_cut = venn_extended(10620000:12620000);
% 
% % Now get the phase of each event type
% event_phase = [];
% index = 1
% while index < length(venn_cut)
%     if venn_cut(index) ~= 0
%         index2 = index;
%         while venn_cut(index2+1) == venn_cut(index2) & index2 < length(venn_cut)
%             index2 = index2 + 1;
%         end
%         mean_phase = mean(phase_degrees(index:index2));
%         event_phase = [event_phase; venn_cut(index) index2-index+1 mean_phase];
%         index = index2+1;
%     else
%         index = index + 1;
%     end
% end
% 
% % Now histogram of degrees per group
% 
% hist_swr = event_phase(event_phase(:,1)==3,3);
% hist_sw = event_phase(event_phase(:,1)==1,3);
% hist_rip = event_phase(event_phase(:,1)==2,3);
% 
% figure
% tiledlayout(1,3)
% t1 = nexttile;
% polarhistogram(deg2rad(hist_swr), 18, 'FaceColor','g');
% title('SWR');
% t2 = nexttile;
% polarhistogram(deg2rad(hist_sw), 18, 'FaceColor','b');
% title('SW');
% t3 = nexttile;
% polarhistogram(deg2rad(hist_rip), 18, 'FaceColor','r');
% title('Ripples');
% 
% %% Analysis during UP/DOWNs
% 
% % For now we restrict to the 2 big NREM periods where UP-DOWNS oscillations
% % settle nicely. Boundaries are:
% bounds = zeros(1,L);
% bounds(4046000:5125000) = 1;
% bounds(10620000:12620000)= 1;
% bounds = bounds == 1;
% 
% PFC203_cut = PFC203(:,1);
% UP_DOWN_cut = UP_DOWN_log;
% PFC203_cut(~bounds) = 0;
% UP_DOWN_cut(~bounds) = 0;
% venn_cut = [repelem(venn,120) zeros(1,L-length(venn)*120)];
% venn_cut(~bounds) = 0;
% venn_cut2 = venn_cut(1:120:end);
% UP_DOWN_cut2 = UP_DOWN_cut(1:120:end);
% 
% figure
% tiledlayout(2,1)
% t1 = nexttile;
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC203_cut)
% for ii=1:length(venn_cut2)
%     if venn_cut2(ii) ~=0
%         patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(PFC203_cut)) -max(abs(PFC203_cut)) max(abs(PFC203_cut)) max(abs(PFC203_cut))], colorsvenn(venn_cut2(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
%         hold on
%     end
% end
% t2 = nexttile;
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L),UP_DOWN_cut)
% for ii=1:length(venn_cut2)
%     if venn_cut2(ii) ~=0
%         patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [0 0 3 3], colorsvenn(venn_cut2(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
%         hold on
%     end
% end
% ylim([0 3])
% linkaxes([t1 t2], 'x');
% 
% % Events during UP or DOWN
% 
% %SWR - venn = 3
% swr_up_down = UP_DOWN_cut2(venn_cut2 == 3);
% 
% figure
% b = bar([0.5 2.5], [sum(swr_up_down==1) sum(swr_up_down==2)]);
% b.FaceColor = 'flat';
% b.EdgeColor = 'none';
% xticks([0.5 2.5])
% xticklabels({'SWR during UP','SWR during DOWN'})
% xtickangle(25)
% b.CData(1,:) = [0.4667    0.7804    0.3529];
% b.CData(2,:) = [0.3412    0.4902    0.1412];
% title('SWR count during UP/DOWN states')
% 
% %RIP alone - venn = 2
% rip_up_down = UP_DOWN_cut2(venn_cut2 == 2);
% 
% figure
% b = bar([0.5 2.5], [sum(rip_up_down==1) sum(rip_up_down==2)]);
% b.FaceColor = 'flat';
% b.EdgeColor = 'none';
% xticks([0.5 2.5])
% xticklabels({'RIP during UP','RIP during DOWN'})
% xtickangle(25)
% b.CData(1,:) = [0.9216    0.2863    0.2863];
% b.CData(2,:) = [0.6000    0.1843    0.1843];
% title('Ripple count during UP/DOWN states')
% 
% %SW alone - venn = 1
% sw_up_down = UP_DOWN_cut2(venn_cut2 == 1);
% 
% figure
% b = bar([0.5 2.5], [sum(sw_up_down==1) sum(sw_up_down==2)]);
% b.FaceColor = 'flat';
% b.EdgeColor = 'none';
% xticks([0.5 2.5])
% xticklabels({'SW during UP','SW during DOWN'})
% xtickangle(25)
% b.CData(1,:) = [0.1255    0.5922    0.9020];
% b.CData(2,:) = [ 0    0.4471    0.7412];
% title('SW count during UP/DOWN states')
% 
% % Look at duration of DOWN states that have an event in them
% mat_duration = [];
% index = 1
% while index < length(UP_DOWN_cut2)
%     if UP_DOWN_cut2(index) == 1
%         index2 = index;
%         while UP_DOWN_cut2(index2) == 1 & index2 < length(UP_DOWN_cut2)
%             index2 = index2 +1;
%         end
%         dura = index2 - index;
%         swr_test = ismember(3, venn_cut2(index:(index2-1)));
%         mat_duration = [mat_duration; dura swr_test];
%         index = index2; 
%     else
%         index = index +1;
%     end
% end
% 
% figure
% h1 = histogram(mat_duration(mat_duration(:,2) == 1, 1)','BinWidth', 1 ,'FaceAlpha', 0.4, 'FaceColor',[0.3922    0.8314    0.0745]);
% hold on
% h2 = histogram(mat_duration(mat_duration(:,2) == 0, 1)', 'BinWidth', h1.BinWidth, 'FaceAlpha', 0.6, 'FaceColor',[0.6510    0.6510    0.6510]);
% hold on
% xline(mean(mat_duration(mat_duration(:,2) == 1, 1)),'--', 'Color',[0.3922    0.8314    0.0745], 'LineWidth', 5)
% hold on
% xline(mean(mat_duration(mat_duration(:,2) == 0, 1)),'--', 'Color',[0.6510    0.6510    0.6510], 'LineWidth', 5)
% title('Distribution of the duration of DOWN states')
% xlabel('Duration (in epoch of 120 ms)')
% legend({'With SWR', 'Without SWR'})
% 
% 
% figure
% boxplot([mat_duration(mat_duration(:,2) == 1, 1)' mat_duration(mat_duration(:,2) == 0, 1)'], [ones(1, length(mat_duration(mat_duration(:,2) == 1, 1))) zeros(1, length(mat_duration(mat_duration(:,2) == 0, 1)))]);
% title('Distribution of the duration of DOWN states')
% xticks([1 2])
% xticklabels({'Without SWR', 'With SWR'})
% a = get(get(gca,'children'),'children');
% set(a(3), 'Color', [0.1490    0.3216    0.0275]); 
% set(a(4), 'Color', [0.3804    0.3804    0.3804]); 
% set(a(5), 'Color', [0.3922    0.8314    0.0745]); 
% set(a(6), 'Color', [0.6510    0.6510    0.6510]); 
% hold on
% yline(mean(mat_duration(mat_duration(:,2) == 1, 1)),'--', 'Color',[0.3922    0.8314    0.0745], 'LineWidth', 2)
% hold on
% yline(mean(mat_duration(mat_duration(:,2) == 0, 1)),'--', 'Color',[0.6510    0.6510    0.6510], 'LineWidth', 2)
% 
% % Phase diagram
% 
% [c,d] = butter(3, [0.1/300 1/300]);
% periodic = filtfilt(c,d,PFC203(:,1));
% 
% figure
% plot(PFC203(:,1))
% hold on
% plot(periodic)
% 
% 
% %% Look at cortical fast oscillations during events
% [c,d] = butter(3, [100/300 299/300]);
% fast_oscillations = filtfilt(c,d,PFC203(:,1));
% 
% %Problem around the removed artifacts ==> set to zero around these regions
% penumbra = 1*600; %ms converted into datapoints, region around artifacts that we remove additionally
% index = 1;
% while index < L
%     if outliers(index) == 1
%         index2 = index;
%         while outliers(index2) == 1
%             index2 = index2+1;
%         end
%         fast_oscillations((index-penumbra):(index2+penumbra-1)) = 0;
%         index = index2+penumbra;
%     else
%         index = index+1;
%     end
% end
% 
% 
% %Get envelope
% env = abs(hilbert(fast_oscillations));
% 
% %Cut to only long NREM
% fast_osci_cut = fast_oscillations;
% fast_osci_cut(PFC203_cut==0) = 0;
% env_cut = env;
% env_cut(PFC203_cut==0) = 0;
% 
% %Detect above threshold
% thresh = 2*std(fast_oscillations(PFC203_cut~=0));
% index_events = env_cut>=thresh;
% 
% %Get start, end and mid of each event
% EVENTS = [];
% index = 1;
% while index < L;
%     if index_events(index) == 1;
%         index2 = index;
%         while index2 < L & index_events(index2) ==1
%             index2 = index2+1;
%         end
%         EVENTS = [EVENTS; index floor(0.5*index + 0.5*(index+2-1)) index2-1];
%         index = index2;
%     else
%         index = index+1;
%     end
% end
% 
% %Postprocessing
% 
% %First, if 2 consecutives events happen within 20ms, consider it's the same
% %event
% gap = 600*20/1000;
% index = 1;
% while index < length(EVENTS)
%     if (EVENTS(index+1,1) - EVENTS(index,3))  <= gap
%         EVENTS(index+1,1) = EVENTS(index,1);
%         EVENTS(index+1,2) = floor((EVENTS(index+1,1) + EVENTS(index+1,3))/2);
%         EVENTS(index,:) = [0 0 0];
%     end
%     index = index+1;
% end
% EVENTS((EVENTS(:,1) == 0),:) = [];
% 
% % Then we remove events that are very short (less than 20 ms)
% outliers = zeros(1,length(EVENTS));
% for ii = 1:length(EVENTS)
%     if (EVENTS(ii,3)-EVENTS(ii,2)) <= gap
%         outliers(ii) = 1;
%     end
% end
% EVENTS(outliers==1,:) = [];
% 
% events = fast_osci_cut;
% events(index == 0) = NaN;
% events(index == 1) = max(fast_osci_cut);
% 
% index_events2 = zeros(L,1);
% for ii = 1:length(EVENTS)
%     index_events2(EVENTS(ii,1):EVENTS(ii,3)) = 1;
% end
% % 
% % figure
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), fast_osci_cut)
% % hold on
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), env_cut)
% % hold on
% % yline(mean(fast_osci_cut)+2*std(fast_oscillations(PFC203_cut~=0)))
% % hold on
% % % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), index_events*20)
% % % hold on
% % plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), index_events2*20)
% 
% 
% % Now we have our events, we want to look at the mean frequency of each
% % event
% 
% frequency_fast_osci = [];
% for ii=1:length(EVENTS)
%     event = fast_osci_cut(EVENTS(ii,1):EVENTS(ii,3));
%     [pxx,f] = pmtm(event,4,[],600);
%     FREQ = meanfreq(pxx, 600, [100 297]);
%     frequency_fast_osci = [frequency_fast_osci; FREQ];
% end
% 
% %Now plot distribution of these frequencies
% 
% figure
% histogram(frequency_fast_osci, 'BinWidth', 2);
% title('Distribution of the mean frequency of fast oscillatory events in the PFC')
% 
% %NB ==> ONLY DURING NREM, because harder to detect events in REM, or short
% %NREM
% %% Multiplets 
% % We want here to detect ripples that are close to one another
% [M_multiplets, Mx_multiplets]=findmultiplets(M);
% timestamps = [M_multiplets.singlets{1} M_multiplets.doublets{1} M_multiplets.triplets{1} M_multiplets.quatruplets{1}];
% % find timestamps that occur more than once, and when they occur, that way
% % we can know which groups are mixed up
% 
% multip = zeros(size(ripplespyr));
% for k = 1:length(timestamps)
%     multip(floor(timestamps(k))*600) = 1*(k<=1104) + 2*(k>1104 & k<=1480) + 3*(k>1480 & k <=1588) + 4 * (k>1588);
% end
% 
% figure
% tiledlayout(2,1)
% t1 = nexttile();
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), ripplespyr)
% 
% t2 = nexttile();
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), multip)
% ylim([-1 5])
% 
% linkaxes([t1 t2], "x")
% 
% %% SWR PCA
% 
% % WAS JUST A FIRST TRY, DIDN'T GET TO INVESTIGATE DEEPER IN THAT DIRECTION
% %We gather by event
% 
% 
% %Now we gather together the non-zero successive elements in venn
% event_index = []; %Will contain the index of beginning/end of each event
% index = 1;
% while index<=tmp
%     if venn(index)~=0
%         start = index;
%         index = index+1;
%         while venn(index)~=0 & index <=tmp
%             index = index+1;
%         end
%         event_index = [event_index [start; index-1]];
%     else
%         index = index+1;
%     end
% end
% 
% 
% %We only run the PCA every "event"
% epdata = data2ep(HPC203(:,1), 0.2, 600);
% [pxx,f] = pmtm(epdata,4,[],400);
% 
% 
% featuresPCA = zeros(12,length(event_index));
% 
% for kk=1:length(event_index)
%     featuresPCA(1,kk) =  peak2peak(HPC203((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)),3)); %For sharp waves, the larger the sharp wave, the larger the value
%     [pxx,f] = pmtm(HPC203((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)),3),4,[],600);
%     featuresPCA(2, kk) = log(spectral_power(90, 200, 1, pxx, f)); %Ripple power
%     featuresPCA(3, kk) = sum(sum(bin_ripples(:,event_index(1,kk):event_index(2,kk)))); %Ripple duration
%     tmp = bin_filt_rip(:,event_index(1,kk):event_index(2,kk));
%     featuresPCA(4,kk) = peak2peak(tmp(:)); %Amplitude of the ripple
%     featuresPCA(5, kk) = log(spectral_power(0.5, 1.5, 1, pxx, f)); 
%     featuresPCA(6, kk) = log(spectral_power(2, 4, 1, pxx, f)); 
%     featuresPCA(7, kk) = log(spectral_power(4, 8, 1, pxx, f)); 
%     featuresPCA(8, kk) = log(spectral_power(8, 20, 1, pxx, f)); 
%     featuresPCA(9, kk) = log(spectral_power(30, 50, 1, pxx, f)); 
%     featuresPCA(10, kk) = log(spectral_power(50, 80, 1, pxx, f)); 
%     featuresPCA(11,kk) = double(sum(ripple_hpc203((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk))))>0)*meanfreq(pxx, 600, [90 200]);% Mean frequency of the ripple in the event (if any, 0 otherwise)
%     featuresPCA(12,kk) = sum(tmpsw((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk)))); %Number of sharp waves in the event
% end
% 
% 
% featuresPCA_norm = featuresPCA ./ max(featuresPCA')';
% 
% [coeff,~,latent,~,explained] = pca(featuresPCA_norm');
% PC1 = sum(featuresPCA_norm .* coeff(:,1), 1);
% PC2 = sum(featuresPCA_norm .* coeff(:,2), 1);   
% 
% PC1 = PC1 ./ max(abs(PC1));
% PC2 = PC2 ./ max(abs(PC2));
% 
% [clusters, C, sumD] = kmeans([PC1' PC2'], 3);
% 
% clus_tot = zeros(1,L);
% for kk=1:length(event_index)
%     clus_tot((1+wl*event_index(1,kk)-wl):(wl*event_index(2,kk))) = clusters(kk);
% end
% clus_tot = clus_tot+1;
% figure
% % scatter(featuresPCA_norm(1,:), featuresPCA_norm(2,:))
% % hold on
% colors = zeros(length(clusters),3);
% colors(clusters == 1,:) = repmat([0 0 1], length(colors(clusters == 1)),1);
% colors(clusters == 2,:) = repmat([1 0 0], length(colors(clusters == 2)),1);
% colors(clusters == 3,:) = repmat([0 1 0], length(colors(clusters == 3)),1);
% colors(clusters == 4,:) = repmat([0 0 0], length(colors(clusters == 4)),1);
% % colors(clusters == 5,:) = repmat([1 0 1], length(colors(clusters == 5)),1);
% scatter(PC1, PC2, [], colors)
% 
% 
%     figure
%     tiledlayout(3,1)
% 
%     t1 = nexttile;
%     N = wl;
%     colors = ['w','b','r','g', 'k', 'm'];
%     index1 = 1;
%     index2 = 1;
%     while index2<L
%         if clus_tot(index2)~= clus_tot(index2+1) || (index2 + 1) == L
%                 cl = clus_tot(index2);
%                 patch([index1 index2 index2 index1], [-max(abs(HPC203(:,3))) -max(abs(HPC203(:,3))) max(abs(HPC203(:,3))) max(abs(HPC203(:,3)))], colors(cl), 'FaceAlpha', 0.5, 'LineStyle', 'none')
%                 hold on
%                 index1 = index2 +1;
%         end
%             index2 = index2 + 1;
% 
%     end
%     plot(HPC203(:,3))
%     title('Clustered ripples ')
%     ylabel('Amplitude')
%     xlabel('Time')
%     
%     t2 = nexttile;
%     index1 = 1;
%     index2 = 1;
%     while index2<L
%         if clus_tot(index2)~= clus_tot(index2+1) || (index2 + 1) == L
%                 cl = clus_tot(index2);
%                 patch([index1 index2 index2 index1], [-max(abs(ripplespyr)) -max(abs(ripplespyr)) max(abs(ripplespyr)) max(abs(ripplespyr))], colors(cl), 'FaceAlpha', 0.2, 'LineStyle', 'none')
%                 hold on
%                 index1 = index2 +1;
%         end
%             index2 = index2 + 1;
% 
%     end
%     plot(ripplespyr)
%     ylabel('Amplitude')
%     xlabel('Time')
%     xlim([0 length(ripplespyr)])
% 
% 
% % Count of each cluster
% count_clusters = [];
% c1 = double(clus_tot==2);
% c2 = double(clus_tot==3);
% c3 = double(clus_tot==4);
% c4 = double(clus_tot==5);
% count_clusters(:,1) = sum(reshape(c1(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600)))); %separate in 1min bins
% count_clusters(:,2) = sum(reshape(c2(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))));
% count_clusters(:,3) = sum(reshape(c3(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))));
% count_clusters(:,4) = sum(reshape(c4(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))));
% 
% tmp = repelem(count_clusters,600*60,1);
% count_clusters = [tmp; zeros(length(ripplespyr) - length(tmp), 4)];
% 
% t3 = nexttile;
% 
% plot(count_clusters(:,1)/wl, 'b')
% hold on
% plot(count_clusters(:,2)/wl, 'r')
% hold on
% plot(count_clusters(:,3)/wl, 'g')
% hold on
% plot(count_clusters(:,4)/wl, 'k')
% title('Count of each cluster')
% 
%     linkaxes([t1 t2 t3], 'x')
%     
%     %% Other rat 209
%     
%     
% patIDhpc2091 = 'F:\rat\cannabis\acutes batch 2\209\2020-04-16_11-29-09\102_CH57.continuous'; %pyramidal layer
% patIDhpc2092 = 'F:\rat\cannabis\acutes batch 2\209\2020-04-16_11-29-09\102_CH63.continuous'; %above pyramidal layer
% patIDhpc2093 = 'F:\rat\cannabis\acutes batch 2\209\2020-04-16_11-29-09\102_CH56.continuous'; %below pyramidal layer
% patIDpfc2091 = 'F:\rat\cannabis\acutes batch 2\209\2020-04-16_11-29-09\102_CH7.continuous'; %channel 5
% patIDpfc2092 = 'F:\rat\cannabis\acutes batch 2\209\2020-04-16_11-29-09\102_CH19.continuous'; %channel 5 from end
% 
% acq_fhz = 30e3;
% ds_fhz = 600;
%  
% Wn = [ds_fhz/acq_fhz ]; % Cutoff=fs_new/2 Hz. 
% [b,a] = butter(3,Wn); %Filter coefficients for LPF.
% 
% HPC209 = [];
% PFC209 = [];
% pathsHPC209 = {patIDhpc2091,patIDhpc2092,patIDhpc2093};
% for ii=1:3
%     [Data, ~, ~] = load_open_ephys_data_faster(pathsHPC209{ii});
%     Data = filtfilt(b,a,Data);
%     HPC209(:,ii) = downsample(Data,acq_fhz/ds_fhz);
% end
% 
% pathsPFC209 = {patIDpfc2091,patIDpfc2092};
% for ii=1:2
%     [Data, ~, ~] = load_open_ephys_data_faster(pathsPFC209{ii});
%     Data = filtfilt(b,a,Data);
%     PFC209(:,ii) = downsample(Data,acq_fhz/ds_fhz);
% end
% 
% TOT = abs(HPC209(:,2)) + abs(HPC209(:,1)) + abs(HPC209(:,3)) + abs(PFC209(:,1))+ abs(PFC209(:,2));
% 
% L = length(HPC209(:,1));
% 
% figure
% plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), TOT)
% 
% tr = 3850; %Visual threshold
% outliers = false(L,1);
% index = 1;
% while index<L
%     if TOT(index)>=tr
%         outliers((index-300):index+1999) = ones(2300,1);
%         index = index+2000;
%     else
%         index = index +1;
%     end
% end
% 
% %Filter out artifacts
% 
% HPC209(outliers,:) = mean(median(HPC209));
% PFC209(outliers,:) = mean(median(PFC209));
% 
% 
% figure
% tiledlayout(5,1)
% tt1 = nexttile;
% plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC209(:,2))
% title('Above pyramidal layer')
% tt2 = nexttile;
% plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC209(:,1))
% title('Pyramidal layer')
% tt3 = nexttile;
% plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC209(:,3))
% title('Below pyramidal layer')
% tt4 = nexttile;
% plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), PFC209(:,1))
% title('Superficial PFC')
% tt5 = nexttile;
% plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), PFC209(:,2))
% title('Deep PFC')
% 
% linkaxes([tt1 tt2 tt3 tt4 tt5], 'x', 'y')