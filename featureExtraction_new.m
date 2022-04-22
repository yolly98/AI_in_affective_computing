%
clc;
clear;
close all;

%%
% Lettura dati
opts = detectImportOptions('evaluations.csv');
opts = setvartype(opts,{'Var1'}, 'string');
Evaluations = readtable('evaluations.csv', opts);
Evaluations.Properties.VariableNames = {'id', 'gender', 'age', 'timestamp', 'emotion', 'level', 'files'};

Dates = table2array(Evaluations(:, 4));
People = table2array(Evaluations(:, 1));
Files = table2array(Evaluations(:, end));

[datesLength, ~] = size(Dates);
% Numero di secondi
windowSize = 10;
matrix = [];
Responses = [];


for p = 1 : datesLength
    disp(p);
    if strcmp(Files{p}, 'NULL') == 0 && length(rmfield(jsondecode(Files{p}), 'time')) > 2
        
        % Lettura file dei sensori
        Muse = readtable(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Muse_0055DAB90EEB/', num2str(Dates(p,:)),'.csv'), 'PreserveVariableNames', true);
        
        if exist(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Shimmer_000666809BE7/', num2str(Dates(p,:)),'.csv'), 'file')
            ShimmerEMG = readtable(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Shimmer_000666809BE7/', num2str(Dates(p,:)),'.csv'), 'PreserveVariableNames', true);
        elseif exist(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Shimmer_000666F38F5E/', num2str(Dates(p,:)),'.csv'), 'file')
            ShimmerEMG = readtable(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Shimmer_000666F38F5E/', num2str(Dates(p,:)),'.csv'), 'PreserveVariableNames', true);
        end
        
        if exist(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Shimmer_000666808EDD/', num2str(Dates(p,:)),'.csv'), 'file')
            ShimmerGSR_PPG = readtable(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Shimmer_000666808EDD/', num2str(Dates(p,:)),'.csv'), 'PreserveVariableNames', true);
        elseif exist(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Shimmer_000666808E2C/', num2str(Dates(p,:)),'.csv'), 'file')
            ShimmerGSR_PPG = readtable(strcat('../SensorsCapture/Person_', num2str(People(p, :)),'/Shimmer_000666808E2C/', num2str(Dates(p,:)),'.csv'), 'PreserveVariableNames', true);
        end
        
        % Determino il minor numero di righe tra i tre file per tener conto di eventuali errori di sincronizzazione
        [museRows, ~] = size(Muse);
        [emgRows, ~] = size(ShimmerEMG);
        [gsrppgRows, ~] = size(ShimmerGSR_PPG);
        minRows = min([museRows, emgRows, gsrppgRows]);
        windowsNumber = fix(minRows/(256*windowSize)) - 1;
        
        for j = 22 : 25
            Muse(:, j) = fillmissing(Muse(:, j), 'linear');
            %Muse(:, j) = filloutliers(Muse(:, j), 'linear');
        end
        
        for i = 1 : windowsNumber
            row = [];
            % Estrazione delle features dai dati EEG
            for j = 22 : 25
                EEGRAW = table2array(Muse(1+i*windowSize*256:(i+1)*windowSize*256, j));
                % Estrazione onde delta, theta, alpha, beta, gamma
                fftEEGRAW = fft(EEGRAW);
                deltaWave = abs(ifft(bandpass(fftEEGRAW, [1 4], 256)));
                thetaWave = abs(ifft(bandpass(fftEEGRAW, [4 7.5], 256)));
                alphaWave = abs(ifft(bandpass(fftEEGRAW, [7.5 13], 256)));
                betaWave = abs(ifft(bandpass(fftEEGRAW, [13 30], 256)));
                gammaWave = abs(ifft(bandpass(fftEEGRAW, [30 44], 256)));
                % Calcolo features
                EEGFeatures = [median(EEGRAW), max(EEGRAW), min(EEGRAW), harmmean(EEGRAW), trimmean(EEGRAW, 10), kurtosis(EEGRAW), skewness(EEGRAW), mean(EEGRAW, 'omitnan'), std(EEGRAW), var(EEGRAW), peak2peak(EEGRAW), peak2rms(EEGRAW), rms(EEGRAW), rssq(EEGRAW), meanfreq(EEGRAW), medfreq(EEGRAW), obw(EEGRAW), max(deltaWave), min(deltaWave), median(deltaWave), mean(deltaWave), max(alphaWave), min(alphaWave), median(alphaWave), mean(alphaWave), max(betaWave), min(betaWave), median(betaWave), mean(betaWave), max(gammaWave), min(gammaWave), median(gammaWave), mean(gammaWave),max(thetaWave), min(thetaWave), median(thetaWave), mean(thetaWave)];
                row = [row, EEGFeatures];
            end
            
            % Estrazione delle features dai dati EMG
            col = find(strcmp(ShimmerEMG.Properties.VariableNames, 'EMG_CH1_24BIT_CAL_mV'), 1);
            ShimmerEMG(:, col) = fillmissing(ShimmerEMG(:, col), 'linear');
			%ShimmerEMG(:, col) = filloutliers(ShimmerEMG(:, col), 'linear');
            
            EMG1 = table2array(ShimmerEMG(1+i*256*windowSize:(i+1)*windowSize*256, col));
            
            col = find(strcmp(ShimmerEMG.Properties.VariableNames, 'EMG_CH2_24BIT_CAL_mV'), 1);
            ShimmerEMG(:, col) = fillmissing(ShimmerEMG(:, col), 'linear');
			%ShimmerEMG(:, col) = filloutliers(ShimmerEMG(:, col), 'linear');
            
            EMG2 = table2array(ShimmerEMG(1+i*256*windowSize:(i+1)*windowSize*256, col));
            
            EMG = EMG2 - EMG1;
            EMGFeatures = [median(EMG), max(EMG), min(EMG), harmmean(EMG), trimmean(EMG, 10), kurtosis(EMG), skewness(EMG), mean(EMG, 'omitnan'), std(EMG), var(EMG), peak2peak(EMG), peak2rms(EMG), rms(EMG), rssq(EMG), meanfreq(EMG), medfreq(EMG), obw(EMG)];
            row = [row, EMGFeatures];%1, EMGFeatures2];
            
            % Estrazione delle features dai dati GSR
            col = find(strcmp(ShimmerGSR_PPG.Properties.VariableNames, 'GSR_Skin_Resistance_CAL_kOhms'), 1);
            ShimmerGSR_PPG(:, col) = fillmissing(ShimmerGSR_PPG(:, col), 'linear');
			%ShimmerGSR_PPG(:, col) = filloutliers(ShimmerGSR_PPG(:, col), 'linear');
            
            GSR = table2array(ShimmerGSR_PPG(1+i*256*windowSize:(i+1)*windowSize*256, col));
            GSRFeatures = [median(GSR), max(GSR), min(GSR), harmmean(GSR), trimmean(GSR, 10), kurtosis(GSR), skewness(GSR), mean(GSR, 'omitnan'), std(GSR), var(GSR), peak2peak(GSR), peak2rms(GSR), rms(GSR), rssq(GSR), meanfreq(GSR), medfreq(GSR), obw(GSR)];
            row = [row, GSRFeatures];
            
            % Estrazione delle features dai dati PPG
            col = find(strcmp(ShimmerGSR_PPG.Properties.VariableNames, 'PPG_A13_CAL_mV'), 1);
            ShimmerGSR_PPG(:, col) = fillmissing(ShimmerGSR_PPG(:, col), 'linear');
			%ShimmerGSR_PPG(:, col) = filloutliers(ShimmerGSR_PPG(:, col), 'linear');
            
            PPG = table2array(ShimmerGSR_PPG(1+i*256*windowSize:(i+1)*windowSize*256, col));
            PPGFeatures = [median(PPG), max(PPG), min(PPG), harmmean(PPG), trimmean(PPG, 10), kurtosis(PPG), skewness(PPG), mean(PPG, 'omitnan'), std(PPG), var(PPG), peak2peak(PPG), peak2rms(PPG), rms(PPG), rssq(PPG), meanfreq(PPG), medfreq(PPG), obw(PPG)];
            row = [row, PPGFeatures];
            
            % Aggiorno matrice aggiungendo la nuova riga
            matrix = [matrix; row];
            Responses = [Responses; Evaluations(p, 5:6)];
        end
    end
end


% Rimpiazzo i valori NAN con 0
%matrix(isnan(matrix)) = 0;
%dataTable = array2table(matrix);
%%
% Nomi delle colonne
variableNames = {'EEG1median','EEG1max','EEG1min','EEG1harmmean','EEG1trimmean','EEG1kurtosis','EEG1skewness','EEG1mean','EEG1std', 'EEG1var','EEG1peak2peak','EEG1peak2rms','EEG1rms','EEG1rssq','EEG1meanfreq','EEG1medfreq','EEG1obw', 'EEG1deltaMax','EEG1deltaMin', 'EEG1deltaMedian','EEG1deltaMean', 'EEG1alphaMax','EEG1alphaMin', 'EEG1alphaMedian','EEG1alphaMean','EEG1betaMax','EEG1betaMin', 'EEG1betaMedian','EEG1betaMean','EEG1gammaMax', 'EEG1gammaMin', 'EEG1gammaMedian','EEG1gammaMean','EEG1thetaMax','EEG1thetaMin', 'EEG1thetaMedian','EEG1thetaMean','EEG2median','EEG2max','EEG2min','EEG2harmmean', 'EEG2trimmean','EEG2kurtosis','EEG2skewness','EEG2mean','EEG2std','EEG2var', 'EEG2peak2peak','EEG2peak2rms','EEG2rms','EEG2rssq', 'EEG2meanfreq','EEG2medfreq','EEG2obw','EEG2deltaMax','EEG2deltaMin', 'EEG2deltaMedian','EEG2deltaMean','EEG2alphaMax','EEG2alphaMin', 'EEG2alphaMedian','EEG2alphaMean','EEG2betaMax','EEG2betaMin', 'EEG2betaMedian','EEG2betaMean','EEG2gammaMax','EEG2gammaMin', 'EEG2gammaMedian','EEG2gammaMean','EEG2thetaMax', 'EEG2thetaMin', 'EEG2thetaMedian','EEG2thetaMean','EEG3median','EEG3max','EEG3min','EEG3harmmean','EEG3trimmean','EEG3kurtosis','EEG3skewness', 'EEG3mean','EEG3std','EEG3var', 'EEG3peak2peak','EEG3peak2rms','EEG3rms','EEG3rssq','EEG3meanfreq','EEG3medfreq','EEG3obw', 'EEG3deltaMax','EEG3deltaMin', 'EEG3deltaMedian','EEG3deltaMean','EEG3alphaMax','EEG3alphaMin', 'EEG3alphaMedian','EEG3alphaMean','EEG3betaMax','EEG3betaMin', 'EEG3betaMedian','EEG3betaMean','EEG3gammaMax','EEG3gammaMin', 'EEG3gammaMedian','EEG3gammaMean','EEG3thetaMax','EEG3thetaMin', 'EEG3thetaMedian','EEG3thetaMean','EEG4median', 'EEG4max','EEG4min','EEG4harmmean','EEG4trimmean','EEG4kurtosis','EEG4skewness','EEG4mean','EEG4std','EEG4var', 'EEG4peak2peak', 'EEG4peak2rms','EEG4rms','EEG4rssq','EEG4meanfreq','EEG4medfreq','EEG4obw','EEG4deltaMax','EEG4deltaMin', 'EEG4deltaMedian','EEG4deltMean', 'EEG4alphaMax','EEG4alphaMin', 'EEG4alphaMedian','EEG4alphaMean','EEG4betaMax','EEG4betaMin', 'EEG4betaMedian','EEG4betaMean','EEG4gammaMax','EEG4gammaMin', 'EEG4gammaMedian','EEG4gammaMean','EEG4thetaMax','EEG4thetaMin', 'EEG4thetaMedian','EEG4thetaMean','EMGmedian','EMGmax','EMGmin','EMGharmmean', 'EMGtrimmean','EMGkurtosis','EMGskewness','EMGmean','EMGstd','EMGvar', 'EMGpeak2peak','EMGpeak2rms','EMGrms', 'EMGrssq','EMGmeanfreq','EMGmedfreq','EMGobw','GSRmedian','GSRmax','GSRmin','GSRharmmean','GSRtrimmean','GSRkurtosis','GSRskewness','GSRmean','GSRstd','GSRvar', 'GSRpeak2peak','GSRpeak2rms','GSRrms','GSRrssq','GSRmeanfreq','GSRmedfreq','GSRobw', 'PPGmedian','PPGmax','PPGmin','PPGharmmean', 'PPGtrimmean','PPGkurtosis','PPGskewness','PPGmean','PPGstd','PPGvar','PPGpeak2peak','PPGpeak2rms','PPGrms','PPGrssq','PPGmeanfreq', 'PPGmedfreq','PPGobw'};
dataTable = array2table(matrix);
dataTable.Properties.VariableNames = variableNames;
Responses.Properties.VariableNames = {'Emotion', 'Level'};

% cols = T.Properties.VariableNames;
% fprintf('Dataset cols with Nan = ');
% colsWithMissing = cols(any(isnan(table2array(dataTable))));
% for i = 1 : size(colsWithMissing, 2)
%     fprintf('%s ', cell2mat(colsWithMissing(i)));
% end
% fprintf('\n');
% col_index = [];
% for i = 1 : size(colsWithMissing, 2)
%     col_index(i) = find(strcmp(cols, colsWithMissing{i}));
% end
% dataTable(:, colsWithMissing) = [];
% cols(col_index) = [];
% 
% outliers = isoutlier(dataTable, 'mean');
% [outliers_row, ~] = find(outliers == 1);
% 
% dataTable(unique(outliers_row), :) = [];
% Responses(unique(outliers_row), :) = [];

T = [dataTable Responses];

% Write dataset to file
if ~exist('../dataset/', 'dir')
    mkdir('../dataset/');
end

writetable(T, '../dataset/data_.csv');


%%
% Tabella normalizzata
dataTableNormalized = normalize(dataTable, 'Range');
TNormalized = [dataTableNormalized, Responses];
% Downsampling e augmentation
if find(strcmp(T.Emotion,'happiness'))
    happiness = downsample(T(find(strcmp(T.Emotion,'happiness')), :), 3);
else
    happiness=[];
end
if find(strcmp(T.Emotion,'boredom'))
    boredom = downsample(T(find(strcmp(T.Emotion,'boredom')), :), 3);
else
    boredom=[];
end
if find(strcmp(T.Emotion,'sadness'))
    sadness = downsample(T(find(strcmp(T.Emotion,'sadness')), :), 2);
else
    sadness=[];
end
if find(strcmp(T.Emotion,'anxiety'))
    anxiety = T(find(strcmp(T.Emotion,'anxiety')), :);
else 
    anxiety=[];
end
if find(strcmp(T.Emotion,'anger'))
    anger = T(find(strcmp(T.Emotion,'anger')), :);
else
    anger=[];
end
if find(strcmp(T.Emotion,'disgust'))
    disgust = T(find(strcmp(T.Emotion,'disgust')), :);
else
    disgust=[];
end
if find(strcmp(T.Emotion,'fear'))
    fear = T(find(strcmp(T.Emotion,'fear')), :);
else
    fear=[];
end
S = RandStream('mt19937ar','Seed',5489);


%%
if ~isempty(anger)
    angerNoise = array2table(awgn(table2array(anger(:, 1:end-2)),10,0,S));
    angerNoise = [angerNoise, anger(:,end-1:end)];
    angerNoise.Properties.VariableNames = anger.Properties.VariableNames;
    reset(S);
    anger = [anger; angerNoise];
end
if ~isempty(disgust)
    for i = 1:2
        S = RandStream('mt19937ar','Seed',5489);
        disgustNoise = array2table(awgn(table2array(disgust(:, 1:end-2)),10,0,S));
        disgustNoise = [disgustNoise, disgust(:,end-1:end)];
        disgustNoise.Properties.VariableNames = disgust.Properties.VariableNames;
        disgust = [disgust;disgustNoise];
        reset(S);
    end
end
if ~isempty(fear)
    S = RandStream('mt19937ar','Seed',5489);
    fearNoise = array2table(awgn(table2array(fear(:, 1:end-2)),10,0,S));
    fearNoise = [fearNoise, fear(:,end-1:end)];
    fearNoise.Properties.VariableNames = fear.Properties.VariableNames;
    fear = [fear;fearNoise];
    reset(S);
end

% Tabella con Augmentation e Downsampling senza normalizzazione
TAugmented=[];
if ~isempty(happiness)
    TAugmented=[TAugmented;happiness];
end
if ~isempty(boredom)
    TAugmented=[TAugmented;boredom];
end
if ~isempty(sadness)
    TAugmented=[TAugmented;sadness];
end
if ~isempty(anxiety)
    TAugmented=[TAugmented;anxiety];
end
if ~isempty(anger)
    TAugmented=[TAugmented;anger];
end
if ~isempty(disgust)
    TAugmented=[TAugmented;disgust];
end
if ~isempty(fear)
    TAugmented=[TAugmented;fear];
end
%TAugmented = [happiness; boredom; sadness; anxiety; anger; disgust; fear];

%% Write dataset to file
if ~exist('../dataset/', 'dir')
    mkdir('../dataset/');
end
writetable(TAugmented, '../dataset/dataset_raw.csv');


%% Select relevant features
opts = detectImportOptions('../dataset/dataset_raw.csv');
opts.SelectedVariableNames = {'GSRmedian', 'GSRmax', 'GSRrms', 'GSRmean' 'PPGmedian', 'PPGmin', 'PPGmax', 'EEG1median', 'EEG1obw', 'EEG1betaMedian', 'EEG2alphaMean', 'EEG3median', 'EEG4median', 'EEG4thetaMin', ...
    'EMGmedian', 'EMGmax', 'EMGmin', 'EMGtrimmean', 'EMGmean', 'EMGrms', 'EMGmedian', 'Emotion'};
datasetF = readtable('../dataset/dataset_raw.csv', opts);


writetable(datasetF, '../dataset/dataset_6emotions.csv');

%estrazione emozioni
testset = datasetF;
 
load("trainedModel.mat");
X = testset(:, 1:end-1);
T = testset(:, end);
 
Y = trainedModel.predictFcn(X);
figure;
plotconfusion(categorical(table2array(T)), Y);