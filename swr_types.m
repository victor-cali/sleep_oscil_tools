% Tuk's script for ripple detection
%% ||||| SWRS TYPES |||||
​
% Script designed to run the analysis of SWRs
% Load various folders with relevants functions afterward
​
%Fieldtrip
addpath('C:\Users\students\Documents\Tugdual\GitHub\fieldtrip');
​
%CorticoHippocampal
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal');
addpath('C:\Users\students\Documents\Tugdual\GitHub\CorticoHippocampal\Ripple_detection');
​
%CBD
addpath('C:\Users\students\Documents\Tugdual\GitHub\CBD');
​
%ADRITOOLS
addpath('C:\Users\students\Documents\Tugdual\GitHub\ADRITOOLS');
​
 
% Additionnal package
addpath('C:\Users\students\Documents\Tugdual\GitHub\bluewhitered'); %Actually I'm not using this one here, it's for a blue/red colormap
addpath('C:\Users\students\Documents\Tugdual\GitHub\analysis-tools');
​
%Chronux
addpath(genpath('C:\Users\students\Documents\Tugdual\chronux'));
%% Rat 203 - vehicle
​
%Paths
patIDhpc2031 = 'E:\rat\cannabis\acutes-batch-2\203\2020-04-02_12-21-50\102_CH59.continuous'; %pyramidal layer
patIDhpc2032 = 'E:\rat\cannabis\acutes-batch-2\203\2020-04-02_12-21-50\102_CH56.continuous'; %above pyramidal layer
patIDhpc2033 = 'E:\rat\cannabis\acutes-batch-2\203\2020-04-02_12-21-50\102_CH48.continuous'; %below pyramidal layer
patIDpfc2031 = 'E:\rat\cannabis\acutes-batch-2\203\2020-04-02_12-21-50\102_CH7.continuous'; %channel 5 in depth
patIDpfc2032 = 'E:\rat\cannabis\acutes-batch-2\203\2020-04-02_12-21-50\102_CH19.continuous'; %channel 5 from end
​
% Acquisition parameters
acq_fhz = 30e3; %acquisition freq
ds_fhz = 600;   %downsampling freq
​
% Design of low pass filter (we low pass to 300Hz)
Wn = [ds_fhz/acq_fhz ]; % Cutoff=fs_new/2 Hz. 
[b,a] = butter(3,Wn); %Filter coefficients for LPF.
​
% Load, filter and store data 
HPC203 = [];
PFC203 = [];
pathsHPC203 = {patIDhpc2031,patIDhpc2032,patIDhpc2033};
for ii=1:3
    [Data, ~, ~] = load_open_ephys_data_faster(pathsHPC203{ii});
    Data = filtfilt(b,a,Data);
    HPC203(:,ii) = downsample(Data,acq_fhz/ds_fhz);
end
​
pathsPFC203 = {patIDpfc2031,patIDpfc2032};
for ii=1:2
    [Data, ~, ~] = load_open_ephys_data_faster(pathsPFC203{ii});
    Data = filtfilt(b,a,Data);
    PFC203(:,ii) = downsample(Data,acq_fhz/ds_fhz);
end
​
TOT = abs(HPC203(:,2)) + abs(HPC203(:,1)) + abs(HPC203(:,3)) + abs(PFC203(:,1))+ abs(PFC203(:,2)); % Sum of amplitudes ==> To visually assess artifacts, as they will appear in every channel and add up
​
L = length(HPC203(:,1));
​
figure
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), TOT)
​
% Artifacts
tr = 5670; %Visual threshold 
outliers = false(L,1);
index = 1;
while index<L
    if TOT(index)>=tr
        outliers((index-300):index+1999) = ones(2300,1); %When we detect an artifact, remove 300 datapoints prior (build up of artifact) and 2000 datapoints after (3.5 sec)
        index = index+2000;
    else
        index = index +1;
    end
end
​
%Filter out artifacts & replace with the mean of the channels' medians
​
HPC203(outliers,:) = mean(median(HPC203)); 
PFC203(outliers,:) = mean(median(PFC203));
​
% Display all prefiltered signals
figure
tiledlayout(5,1)
tt1 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC203(:,2))
title('HPC - above pyramidal layer')
tt2 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC203(:,1))
title('HPC - pyramidal layer')
tt3 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), HPC203(:,3))
title('HPC - below pyramidal layer')
tt4 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), PFC203(:,1))
title('PFC - shallow layer')
tt5 = nexttile;
plot(linspace(duration([0 0 0]),duration([0 0 L/600]),L), PFC203(:,2))
title('PFC - deep layer')
​
linkaxes([tt1 tt2 tt3 tt4 tt5], 'x', 'y')
​
%% Sharp wave ripples detection
% For the detection of SW (sharp waves) we use the BPL (Below Pyramidal
% Layer) of the HPC
%Band pass filter to better see the sharp waves
[c,d] = butter(3, [1/300 20/300]);
filtbpl = filtfilt(c,d,HPC203(:,3));
spw = double(filtbpl <= mean(filtbpl)-5*std(filtbpl)); % Detect sharp waves as signal lower than ?-5*SD
dspw = abs(diff(spw)); % Use absolute derivative for practical purposes: every sharp wave will hence start with a 1 and end with a 1
​
% Detect sharp waves in the HPC layer below pyramidal layer (BPL)
​
%Get the local peak of each sharp wave 
spwpeak = nan(L,1);
index = 1;
while index < L
    if dspw(index) == 1
        index2 = index;
        index = index+1;
        while dspw(index) == 0
            index = index + 1;
        end
        locmin = min(HPC203((index2+1):(index+1),3));
        indmin = find(HPC203(:,3) == locmin);
        spwpeak(indmin) = locmin;
        index = index+1;
    else
        index = index+1;
    end
end
​
% Sharp Waves characterization
N_SW = sum(spwpeak); %Number of SW
sw_epoch = ~isnan(spwpeak);
sw_epoch = reshape(sw_epoch(1:floor(L/(60*600))*60*600), 60*600,floor(L/(60*600))); %separate in 1min bins
sw_fhz = sum(sw_epoch); %Nb of SW per second for each 1min bin
figure
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),length(sw_fhz)), sw_fhz)
title('Sharp Wave rate (/sec) per 1 min bin')
% What we observe: bursts in SW activity for every NREM episode, with an
% apparent decrease as the bout lasts
​
%Ripple detection 
%For ripple detection, we use the pyramidal layer of the HPC
%First we bandpass for ripple frequency range
[e,f] = butter(3, [90/300 200/300]);
ripplespyr = filtfilt(e,f,HPC203(:,1));
​
yourtimevector = (1:length(ripplespyr))/600;
​
thresh = mean(ripplespyr) + 5*std(ripplespyr); %threshold for ripple detection
[S, E, M] = findRipplesLisa(ripplespyr', yourtimevector, thresh, (thresh)*(1/2), [] ); %Adrian's detection
​
​
M_dur = NaT(1, length(M)) - NaT(1);   % Not really important
for kk=1:length(M)
    M_dur(kk) = duration([0 0 M(kk)]);
end
​
% We have the start/end time of each ripple, we create one long logical
% vector where 1's are our ripples
ripple_hpc203 = zeros(1,L);
​
index = 1; % Timestamp index.
index2 = 1; %Events index.
while index < length(yourtimevector) & index2<length(E)
    if S(index2) == yourtimevector(index)
        while E(index2) ~= yourtimevector(index)
            ripple_hpc203(index) = 1;
            index = index + 1;
        end
        index2 = index2 + 1;
    else
        index = index + 1;
    end
end
​
% Now we create a "Venn diagram" of ripples alone, SW alone and co-occuring
% SWripples (or SWR)
% We group by bins of 200ms
wl = 600*0.2; %200ms bins ==> 600 points = 1sec
tmp = floor(L/wl);
​
bin_PL = reshape(HPC203(1:tmp*wl,1), wl, tmp); %Pyramidal Layer
bin_BPL = reshape(HPC203(1:tmp*wl,3), wl, tmp); % Below Pyramidal Layer
bin_ripples = reshape(ripple_hpc203(1:tmp*wl), wl, tmp); 
bin_filt_rip = reshape(ripplespyr(1:tmp*wl), wl, tmp); % Filtered pyramidal layer for ripples
tmpsw = ~isnan(spwpeak);
bin_sw = reshape(tmpsw(1:tmp*wl), wl, tmp);
​
% Sum per bin ==> if sum is >0, means that one event is happening (either
% SW in the BPL or ripple in the pyramidal layer - PL)
sumbin_sw = sum(bin_sw); 
sumbin_rip = sum(bin_ripples);
​
venn = zeros(1,tmp);
for kk=1:tmp
    venn(kk) = 1*(sumbin_sw(kk)>0 & sumbin_rip(kk) == 0) + 2*(sumbin_sw(kk) == 0 & sumbin_rip(kk) > 0) + 3*(sumbin_sw(kk) > 0 & sumbin_rip(kk) > 0);
end
% So in the end that is: 1 for SWs alone, 2 for ripples alone, 3 for SWR
​
​
COUNT = zeros(1,3);
​
% Need to adjust venn for consecutive "events" = one longer event
index1 = 1;
while index1 < length(venn)
    if venn(index1) ~=0 % Beginning of an "event"
       index2 = index1+1;
       while venn(index2) ~=0 & index2 < length(venn) %We look the length of that event (or consecutive events)
           index2 = index2+1;
       end
       event = venn(index1:(index2 -1));
       if length(event) >= 2 %If that event last longer than 2 200ms-bins, then it potentially is a mix of different "types" of events
           % Several cases
           if ismember(3, event)
               venn(index1:(index2-1)) = 3; %If one 200 ms bin is a SWR, then the whole event is a SWR
           elseif ismember(1, event) & ismember(2, event) 
                venn(index1:(index2-1)) = 3; %If it contains both a RIP without SW and a SW without RIP, then it's a SWR
           end
       end
       COUNT = COUNT + ismember(1,venn(index1:(index2-1)))*[1 0 0] + ismember(2,venn(index1:(index2-1)))*[0 1 0]+ ismember(3,venn(index1:(index2-1)))*[0 0 1];
       index1 = index2;
    else
        index1 = index1 +1;
    end
end
​
% Barplot of the number of events
figure
b = bar([0.5 2 3.5], COUNT);
b.FaceColor = 'flat';
b.EdgeColor = 'none';
xticks([0.5 2 3.5])
xticklabels({'SW without ripple','Ripple without SW', 'SWR'})
xtickangle(25)
b.CData(1,:) = [0    0.4471    0.7412];
b.CData(2,:) = [0.9216    0.2863    0.2863];
b.CData(3,:) = [0.4667    0.7804    0.3529];
title('Number of each event type for rat 203')
​
% Plots
% Note - Patches can be really long to display in plots
figure
tiledlayout(4,1)
% t1 = nexttile;
colorsvenn = ['b', 'r', 'g'];
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), HPC203(:,1))
% hold on
% for ii=1:length(venn)
%     if venn(ii) ~=0
%         patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(HPC203(:,1))) -max(abs(HPC203(:,1))) max(abs(HPC203(:,1))) max(abs(HPC203(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
%         hold on
%     end
% end
% title('Pyramidal Layer')
t1 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), ripplespyr)
hold on
stem(M_dur,600*ones(size(M_dur)))
title('Pyramidal layer filtered for ripples')
t2 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), HPC203(:,3))
hold on
% plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), filtblwpyr)
% hold on
scatter(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L),spwpeak, 'red', 'filled')
title('Below pyramidal layer filtered with sharp waves')
t3 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), PFC203(:,1))
hold on
for ii=1:length(venn)
    if venn(ii) ~=0
        patch([(1+(ii-1)*wl) (ii*wl) (ii*wl) (1+(ii-1)*wl)]*1/(600*24*3600), [-max(abs(HPC203(:,1))) -max(abs(HPC203(:,1))) max(abs(HPC203(:,1))) max(abs(HPC203(:,1)))], colorsvenn(venn(ii)), 'FaceAlpha', 0.5, 'LineStyle', 'none' )
        hold on
    end
end
title('PFC channel')
t4 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS'),L), [0; diff(PFC203(:,1))])
linkaxes([t1 t2 t3 t4], 'x')
​
​
%% Amplitudes of sharp waves vs ripples
% filtbpl for the sharp waves
% ripplespyr for the ripples
amplitudes = [];
for ii=1:length(venn)
   if venn(ii)~=0 
       sw_vec = filtbpl(((ii-1)*120+1):ii*120); % 1 datapoint in venn = 120 datapoints in the raw signal (200ms * 600 datapoints/sec)
       rip_vec = ripplespyr(((ii-1)*120+1):ii*120);
       pfc1 = PFC203(((ii-1)*120+1):ii*120,1);
       pfc2 = PFC203(((ii-1)*120+1):ii*120,2);
       ampl_sw = mean(filtbpl)-min(sw_vec);
       ampl_rip = peak2peak(rip_vec);
       ampl_pfc1 = peak2peak(pfc1);
       ampl_pfc2 = peak2peak(pfc2);
       amplitudes = [amplitudes; [ampl_sw ampl_rip ampl_pfc1 ampl_pfc2]];
   end
end
​
figure
scatter(amplitudes(:,1), amplitudes(:,2))
title('Amplitude of sharp waves vs ripples')
xlabel('Sharp wave amplitude') 
ylabel('Ripple amplitude')
​
%No clear tendency, although there seems to be some kind of correlation
%% REMINDER
% Venn: 1=SW, blue, 2=rip, red, 3=SWR, green
%% Spectral power in PFC after events 
​
% Spectrogram, using an overlap window of 1sec (timestep 0.2 sec)
params=struct('tapers',[],'Fs',[],'fpass',[]);
params.tapers = [5/2 4];
params.Fs = 600;
params.fpass = [1 200];
[S,t,f] = mtspecgramc(PFC203(:,1),[1 0.2],params);
​
figure
tiledlayout(2,1)
t1 = nexttile;
surf(duration(0,0,t),f,10*log10(S'),'EdgeColor', 'none')
view(2)
colormap(t1, 'parula')
%colorbar();
set(gca,'YDir','normal')
title('Spectrogram')
​
t2 = nexttile;
plot(linspace(duration(0,0,0,0, 'Format', 'hh:mm:ss.SSS'),duration(0, 0, L/600, 0, 'Format', 'hh:mm:ss.SSS...