clc
clear all
close all

%% Load ECG and clinical information
ECG = load("Example_ECG_signal.mat");
ECG_signal = double(ECG.pEF_1008);
clinical_information_table = readtable("Clinical_information.csv");
clinical_information = table2array(clinical_information_table(1,2:14));
fs = clinical_information_table.SamplingFreq;
EF = clinical_information_table.EF;

%% Extract HRV
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ECG_signal,fs,0);
RR_time = diff(qrs_i_raw)./fs;
%%%%%% Smooth and remove error peaks
t = 1:length(RR_time);
T= [0.07,0.085,0.15,0.5];
[Xd,Xnoise,NoiseLoc] = SDROM(RR_time,9,T);
RR_time2 = ada_f(Xd);

%%%%%% Find starting hour
time_all = double.empty;
RSS_all = double.empty;
cum_time = cumsum(RR_time2);
data = RR_time2;
time = cum_time;
Kstep = 30;
meanData=[];
meanTime=[];
k_old = 1;
for k=0:Kstep:length(data)-Kstep
    dataWindow=data(k_old:k+Kstep);
    meanData=[meanData mean(dataWindow)];
    k_old = k+Kstep;
    meanTime=[meanTime time(k+floor(Kstep/2))];
end
meanData_seconds = meanData;
cum_time_edited = cumsum(meanData_seconds);
[r c] = size(meanData_seconds);
pq = (0:c-1)./10;
t = interp1(double(cum_time_edited),pq,'linear');     
w = 2*pi;
[M, Amp, phi, RSS] = cosinor(t(11:end)./fs,meanData_seconds(11:end),w,0.05);
degree = (phi*180)/pi;
time = abs((degree/360)*24);
time_rounded = round(time);

%%%%%% Extract per hour
cum_time = cumsum(RR_time2);
cum_time(cum_time > 86400) = [];
time_old = 1;
number = 0;
RR_time_1hour = cell.empty;
    for j = 3600:3600:86400
    cum_time_1hour_more = find(cum_time >= j);
    if isempty(cum_time_1hour_more) == 0
    cum_time_1hour_first = cum_time_1hour_more(1);
    cum_time_1hour{number+1,1} = cum_time(time_old:cum_time_1hour_first);
    RR_time_1hour{number+1,1} = RR_time2(time_old:cum_time_1hour_first);
    time_old = cum_time_1hour_first;
    number = number + 1;
    else
    cum_time_1hour{number+1,1} = cum_time(cum_time_1hour_first:end);
    RR_time_1hour{number+1,1} = RR_time2(cum_time_1hour_first:end);
    end
    end

%%%%%% Adjust hours to start from 00:00
RR_added_hour_cell = cell.empty;
RR_patient_fixed = cell.empty;
RR_time_split_fixed = cell.empty;
for patientID = 1:1
    RR_patient = RR_time_1hour;
    Starting_time = time_rounded;
    added_hours = 24 - round(Starting_time);
    if added_hours > 1
       RR_patient_fixed = RR_patient;
       RR_length = length(RR_patient);
       RR_missing = 24 - RR_length;
       if RR_missing ~= 0
         for RR_ID = 1:RR_missing
       RR_patient_fixed = [RR_patient_fixed; NaN];
       end
    end
    for hourID = 1:added_hours
        RR_added_hour = RR_patient_fixed{hourID,1}; %Per hour
        RR_added_hour_cell{1,1} = RR_added_hour;
        RR_patient_fixed(hourID) = {[]};
        RR_patient_fixed = [RR_patient_fixed;RR_added_hour_cell];
    end
       RR_patient_fixed(cellfun('isempty',RR_patient_fixed)) = [];
       RR_time_split_fixed{patientID,1} = RR_patient_fixed;
       RR_patient = cell.empty;
       RR_added_hour_cell = cell.empty;
       RR_patient_fixed = cell.empty;
    else
    RR_time_split_fixed{patientID,1} = RR_patient;
    end
[sizer sizec] = size(RR_time_split_fixed{patientID,1});
missing_hours = 24-sizer;
if sizer ~= 24
add = RR_time_split_fixed{patientID,1};
    for count = 1:missing_hours
        add = [add;NaN];
    end
RR_time_split_fixed{patientID,1} = add;
end
end
RR_time_1hour_fixed = RR_time_split_fixed{1,1};

%% Extract HRV features
for hour_id = 1:24   
selected_hour = RR_time_1hour{hour_id,1};
[hrv_td] = hrv_time_edit(selected_hour./1000,50/1000);
hrv_td2 = table2array(hrv_td);
time_features = hrv_td2';
                 
[hrv_fd] = hrv_freq_edit(selected_hour);
hrv_fd2 = table2array(hrv_fd);
freq_features = hrv_fd2';
                
[hrv_nl] = hrv_nonlinear_edit(selected_hour./1000);
hrv_nl2 = table2array(hrv_nl);
nonlinear_features = hrv_nl2';
                
[hrv_frag] = hrv_fragmentation_edit(selected_hour./1000);
hrv_frag2 = table2array(hrv_frag);
fragmentation_features = hrv_frag2';
         
all_hrv_features(:,hour_id) = [time_features;freq_features;nonlinear_features;fragmentation_features];
end  
                 
load mean_features %Based on the paper data for normalization
load std_features %Based on the paper data for normalization
                 
selected_feature2 = (all_hrv_features-mean_features)./std_features;
bottom = 1+exp((-selected_feature2));
top = 1;
all_hrv_features_normalized = top./bottom;

%% Generate HRV images
[features_set,new3] = generate_images(clinical_information,all_hrv_features_normalized);

%%%%%% Select best performing features (RMSSD, SEM, HFNORM, HFPOWER, LFNORM, SD1, IALS) 
features_set_best = features_set(:,:,[3;17;5;9;7;10;23]);

figure('Position',[608,345,632,533]);
imshow(features_set_best(:,:,1))
black = [0 0 0]./255;
teal = [22 149 146]./255;
purple = [122 0 244]./255;
dark_orange = [252 133 14]./255;
white = [255 255 255]./255;
mycolormap = customcolormap([0 0.25 0.5 0.75 1], [white;dark_orange;purple;teal;black]);
colormap(mycolormap)
caxis([0 1]);
colorbar('FontSize',15,'Ticks',[0:0.1:1]);

%% Predict heart failure stage using deep learning
HF_stages = {'HFpEF', 'HFmEF', 'HFrEF'};
network = load('network.mat');
network = SeriesNetwork.loadobj(network.network);
[predictedLabel scores] = classify(network,features_set_best);

['Predicted stage: ', HF_stages{predictedLabel}, ', Score: ', num2str(round(scores(predictedLabel),3))]

%%%%%%%%%%%%%%%%%%%%%%%%%%%% NETWORK
%%%% Network (can be used to train your own data)
% layers = [
%     imageInputLayer(imageSize,'Normalization','zerocenter')
% 
%     convolution2dLayer(3,32,'Stride',2,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
% 
%     groupedConvolution2dLayer(3,1,'channel-wise','Stride',2,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
% 
%     convolution2dLayer(1,32,'Stride',2,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
% 
%     maxPooling2dLayer(2,'Stride',2,'Padding','same')
% 
%     fullyConnectedLayer(3)    
%     softmaxLayer
%     classificationLayer('Classes',classes_weights,'ClassWeights',classWeights,'Name','classification')];
% 
% options = trainingOptions('sgdm',...
%     'ExecutionEnvironment','gpu',...
%     'Minibatchsize',128,...
%     'MaxEpochs',45,...
%     'InitialLearnRate',0.001,...
%     'L2Regularization',0.0001,... 
%     'Shuffle','every-epoch',...
%     'Verbose',false,...
%     'Plots','training-progress');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Interpret the prediction
map = occlusionSensitivity(network,features_set_best,categorical(predictedLabel),...
                           'MaskSize','auto','MiniBatchSize',128);
map = rescale(map,0,1);

%%%%%% Plot heatmap
figure('Position',[608,345,632,533]);
imshow(map);
colormap jet;
colorbar('FontSize',18,'Location','eastoutside','Ticks',[0:0.1:1]);
caxis([0 1]);

figure('Position',[608,345,632,533]);
imshow(features_set_best(:,:,1));
hold on;
imagesc(map,'AlphaData',0.35);
colormap jet;
colorbar('FontSize',18,'Location','eastoutside','Ticks',[0:0.1:1]);
caxis([0 1]);

%%%%%% Plot time importance
angles_radian = (0:360/24:360).*(pi/180);
angles_radian = angles_radian-0.52;
rho = 255.*ones(1,25);
[x,y] = pol2cart(angles_radian,rho);
cart = [x;y];
cart(1,:) = round(cart(1,:)+256);
cart(2,:) = abs(round(cart(2,:)-256));
mask = cell.empty;
for i=1:1:24
    if i == 24
       x = [255,cart(1,i)+0.5,cart(1,1)+0.5];
       y = [255,cart(2,i),cart(2,1)];
       mask{i,1} = poly2mask(x, y, 512, 512);
    end
x = [255,cart(1,i)+0.5,cart(1,i+1)+0.5];
y = [255,cart(2,i),cart(2,i+1)];
mask{i,1} = poly2mask(x, y, 512, 512);
end
array = [8,7,6,5,4,3,2,1,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9];
mask_edited = mask(array);

A_total = double.empty;
A_total2 = cell.empty;
max_values = double.empty;
length_all = double.empty;
max_values_location = double.empty;
count = 0;
for q=1:24   
A = mask_edited{q,1}.*map;    
A_fixed = A(:);
A_total2{q,1} = A_fixed;
length_all(q,1) = length(A_fixed);
[maxvalue,maxvalue_id] = max(A_fixed);
A_total = [A_total;A_fixed];
max_values(q,1) = maxvalue;
if q == 1
   max_values_location(q,1) = maxvalue_id+(0*count);
else 
max_values_location(q,1) = maxvalue_id+(length_all(q-1)*count);
end

count = count + 1;
end
[final_max_value_id,final_max_value] = find(max_values >= 0.7);
final_max_values_location = max_values_location(final_max_value_id,1);

figure('Position',[358,249,1243,603]);
plot([1:length(A_total)],A_total,'-k','LineWidth',2)
hold on
ax = gca;
ax.FontSize = 20;
xlabel('Time (h)','FontWeight','bold')
ylabel('Region importance','FontWeight','bold')
xlim([0 length(A_total)])
xticks([0:length(A_total)/24:length(A_total)])
xticklabels({'','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24'})
xtickangle(0)
grid
[up,lo] = envelope(A_total,200000,'peak');
plot([1:length(A_total)],up,'-','linewidth',2.5,'Color',[0.00,0.45,0.74])
plot([length(A_total)*2:(length(A_total)*2)+1],[1 1],'-r','LineWidth',2)
plot([1:length(A_total)],0.7*ones(1,length(A_total)),'--r','LineWidth',2)
xlim([0 length(A_total)])
ylim([0 1.05])
yticks([0:0.1:1])

for o = 1:length(final_max_value_id)
    selected_id = final_max_value_id(o,1);
    start = length(A_total2{1,1})*selected_id - length(A_total2{1,1}) + 1;
    plot([start:length(A_total2{selected_id,1})*selected_id],A_total2{selected_id,1},'-r','LineWidth',2)
end
legend('Heatmap regions','Upper envelope','Important regions','Importance threshold',...
       'FontWeight','normal','orientation','horizontal','location','northoutside','FontSize',17);

final_max_value_id0 = final_max_value_id;
%%%%%% Plot clinical information importance
timeset = final_max_value_id0;

figure('Position',[358,249,1243,603]);
hold on
ax = gca;
ax.FontSize = 20;
xlabel('Demographic and clinical features','FontWeight','bold')
ylabel('Region importance','FontWeight','bold')
xticklabels({'','Age','Sex','BMI','Smoking','Diabetes','HT','AP','VT','Prior MI','Beta-b.','ACE-inhib.','Anti-arr.','Diuretic'})
xtickangle(60)
ylim([0 1.05])
yticks([0:0.1:1])
grid    
box
color_increase = 0;
for i = 1:length(timeset)
    selected_time = timeset(i);
    selected_time_mask = mask_edited{selected_time,1};
    
        masking_heatmap_perhour = cell.empty;
    for j = 1:length(new3)
        selected_circle_mask = logical(new3{j,1});
        selected_time_circle_mask = selected_circle_mask.*selected_time_mask;
        masking_heatmap_perhour{j,1} = selected_time_circle_mask.*map;
    end

A_total_demog = double.empty;
A_total2_demog = cell.empty;
max_values_demog = double.empty;
length_all_demog = double.empty;
max_values_location_demog = double.empty;
count_demog = 0;
for q=1:size(masking_heatmap_perhour,1)  
A_demog = masking_heatmap_perhour{q,1};    
A_fixed_demog = A_demog(:);
A_total2_demog{q,1} = A_fixed_demog;
length_all_demog(q,1) = length(A_fixed_demog);
[maxvalue,maxvalue_id] = max(A_fixed_demog);
A_total_demog = [A_total_demog;A_fixed_demog];
max_values_demog(q,1) = maxvalue;
if q == 1
   max_values_location_demog(q,1) = maxvalue_id+(0*count_demog);
else 
max_values_location_demog(q,1) = maxvalue_id+(length_all_demog(q-1)*count_demog);
end

count_demog = count_demog + 1;
end
[final_max_value_id,final_max_value] = find(max_values_demog >= 0.7);
final_max_values_location_demog = max_values_location_demog(final_max_value_id,1);

[up,lo] = envelope(A_total_demog,200000,'peak');

if double(predictedLabel) == 1
    Color = [0.25 0.36 0]; %HFpEF
end
if double(predictedLabel) == 2
    Color = [0 0 0.28]; %HFmEF
end
if double(predictedLabel) == 3
    Color = [0.3 0 0.11]; %HFrEF
end

Color(2) = Color(2)+color_increase;

plot([1:length(A_total_demog)],up,'-','linewidth',2.5,'Color',Color)

color_increase = color_increase + 0.08;

end
xlim([0 length(A_total_demog)])
xticks([0:length(A_total_demog)/13:length(A_total_demog)])

plot([1:length(A_total_demog)],0.7*ones(1,length(A_total_demog)),'--r','LineWidth',2)

%%%%%%%%%% LEGEND GOES FROM DARK TO LIGHT AS HOUR 00:00 onwards
