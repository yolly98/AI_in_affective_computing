format longG
% Lettura dati
opts = detectImportOptions("evaluations.csv");
opts = setvartype(opts,{'Var1'},'string');
Evaluations = readtable("evaluations.csv",opts);
Dates = table2array(Evaluations(:,4));
People = table2array(Evaluations(:,1));
[datesLength,~] = size(Dates);
windowSize = 10; % Numero di secondi
formatSpec = '%.0f';
matrix = [];
Responses = [];
for p = 1:datesLength %per ogni entrata di evaluations
  disp(p)
  % Lettura file dei sensori
  MusePath=strcat("../SensorsCapture/Person_", num2str(People(p, :),formatSpec),"/Muse_0055DAB90EEB/", string(Dates(p,:)),".csv");
  Shimmer1Path=strcat("../SensorsCapture/Person_", num2str(People(p, :),formatSpec),"/Shimmer_000666809BE7/", string(Dates(p,:)),".csv");
  Shimmer2Path=strcat("../SensorsCapture/Person_", num2str(People(p, :),formatSpec),"/Shimmer_000666808E2C/", string(Dates(p,:)),".csv");
  
  if ~exist(MusePath,"file") 
      disp(strcat("Not found: ",MusePath)); 
      continue;
  end
  Muse = readtable(MusePath,'ReadVariableNames',false);
  
  if ~exist(Shimmer1Path,"file")
      disp(strcat("Not found: ",Shimmer1Path)); 
      continue;
  end
  ShimmerEMG = readtable(Shimmer1Path,'ReadVariableNames',false);
  
  if ~exist(Shimmer2Path,"file") 
      disp(strcat("Not found: ",Shimmer2Path)); 
      continue;
  end
  ShimmerGSR_PPG = readtable(Shimmer2Path,'ReadVariableNames',false);
  
  % Determino il minor numero di righe tra i tre file per tener conto di eventuali errori di sincronizzazione
  [museRows,~] = size(Muse);
  [emgRows,~] = size(ShimmerEMG);
  [gsrppgRows,~] = size(ShimmerGSR_PPG);
  minRows = min([museRows,emgRows,gsrppgRows]);
  windowsNumber = fix(minRows/(256*windowSize))-1;
  
  if(windowsNumber<=0)
      disp(strcat("registrazione ",string(Dates(p,:))," non valida"));
      continue;
  end
  
  for i = 0:windowsNumber  %per ogni misurazione
    disp(strcat("windowsNumber: ",num2str(i)));
    row = [];
    
    % Estrazione delle features dai dati GSR
    disp("estrazione GSR");
    GSR = table2array(ShimmerGSR_PPG(1+i*256*windowSize:(i+1)*windowSize*256,5));
    GSRFeatures = [median(GSR), max(GSR), rms(GSR), mean(GSR, 'omitnan')];
    row=[row,GSRFeatures];
    
    % Estrazione delle features dai dati PPG
    disp("estrazione PPG");
    PPG = table2array(ShimmerGSR_PPG(1+i*256*windowSize:(i+1)*windowSize*256,3));
    PPGFeatures = [median(PPG), min(PPG), max(PPG)];
    row=[row,PPGFeatures];
    
    % Estrazione delle features dai dati EEG
    disp("estrazione EEG");
    for j = 22:25
      %disp(Muse)
      EEGRAW = table2array(Muse(1+i*windowSize*256:(i+1)*windowSize*256,j));
      %disp(EEGRAW)
      % Estrazione onde delta, theta, alpha, beta, gamma
      if(j==22)
          fftEEGRAW = fft(EEGRAW);
          betaWave = abs(ifft(bandpass(fftEEGRAW, [13 30], 256)));
          EEGFeatures = [median(EEGRAW), obw(EEGRAW), median(betaWave)];
      elseif(j==23)
          fftEEGRAW = fft(EEGRAW);
          alphaWave = abs(ifft(bandpass(fftEEGRAW, [7.5 13], 256)));
          EEGFeatures = mean(alphaWave);
      elseif(j==24)
          EEGFeatures = median(EEGRAW);
      elseif(j==25)
          fftEEGRAW = fft(EEGRAW);
          thetaWave = abs(ifft(bandpass(fftEEGRAW, [4 7.5], 256)));
          EEGFeatures = [median(EEGRAW), min(thetaWave)];
      end
      
      row=[row,EEGFeatures];
    end
    
    % Estrazione delle features dai dati EMG
    disp("estrazione EMG");
    EMG1 = table2array(ShimmerEMG(1+i*256*windowSize:(i+1)*windowSize*256,4));
    EMG2 = table2array(ShimmerEMG(1+i*256*windowSize:(i+1)*windowSize*256,5));
    EMG=[];
    for j=1:length(EMG1)
        EMG=[EMG;EMG1(j)-EMG2(j)];
    end
    EMGFeatures = [median(EMG), max(EMG), min(EMG), trimmean(EMG, 10), mean(EMG, 'omitnan'), rms(EMG)];
    row=[row,EMGFeatures];
   
    % Aggiorno matrice aggiungendo la nuova riga
    matrix=[matrix;row];
    Responses = [Responses;Evaluations(p,5:6)];
  end
end
disp("termine estrazione");
% Rimpiazzo i valori NAN con 0
matrix(isnan(matrix))=0;
dataTable = array2table(matrix);
% Nomi delle colonne
variableNames={'GSRmedian','GSRmax','GSRrms','GSRmean','PPGmedian','PPGmin','PPGmax','EEG1median','EEG1obw','EEG1betaMedian','EEG2alphaMean','EEG3median','EEG4median','EEG4thetaMin','EMGmedian','EMGmax','EMGmin','EMGtrimmean','EMGmean','EMGrms'};
dataTable.Properties.VariableNames = variableNames;
Responses.Properties.VariableNames = {'Emotion', 'Level'};
T = [dataTable, Responses];
% Tabella normalizzata
dataTableNormalized = normalize(dataTable, 'Range');
TNormalized = [dataTableNormalized, Responses];

%------------------------- Downsampling e augmentation
if isempty(find(strcmp(T.Emotion,'happiness')))
    happiness=[];
else
    happiness = downsample(T(find(strcmp(T.Emotion,'happiness')), :), 3);
end
if isempty(find(strcmp(T.Emotion,'boredom')))
    boredom=[];
else
    boredom = downsample(T(find(strcmp(T.Emotion,'boredom')), :), 3);
end
if isempty(find(strcmp(T.Emotion,'sadness')))
    sadness=[];
else
    sadness = downsample(T(find(strcmp(T.Emotion,'sadness')), :), 2);
end
if isempty(find(strcmp(T.Emotion,'anxiety')))
    anxiety=[];
else
    anxiety = T(find(strcmp(T.Emotion,'anxiety')), :);
end
if isempty(find(strcmp(T.Emotion,'anger')))
    anger=[];
else
    anger = T(find(strcmp(T.Emotion,'anger')), :);
end
if isempty(find(strcmp(T.Emotion,'disgust')))
    disgust=[];
else
    disgust = T(find(strcmp(T.Emotion,'disgust')), :);
end
if isempty(find(strcmp(T.Emotion,'fear')))
    fear=[];
else
    fear = T(find(strcmp(T.Emotion,'fear')), :);
end
if ~isempty(anger)
    S = RandStream('mt19937ar','Seed',5489);
    angerNoise = array2table(awgn(table2array(anger(:, 1:20)),10,0,S));
    angerNoise = [angerNoise, anger(:,21:22)];
    angerNoise.Properties.VariableNames = anger.Properties.VariableNames;
    reset(S);
else
    angerNoise=[];
end
anger = [anger;angerNoise];
  
if ~isempty(disgust)
    for i = 1:2
      S = RandStream('mt19937ar','Seed',5489);
      disgustNoise = array2table(awgn(table2array(disgust(:, 1:20)),10,0,S));
      disgustNoise = [disgustNoise, disgust(:,21:22)];
      disgustNoise.Properties.VariableNames = disgust.Properties.VariableNames;
      reset(S);
    end
else
    disgustNoise=[];
end
disgust = [disgust;disgustNoise];
if ~isempty(fear)
    S = RandStream('mt19937ar','Seed',5489);
    fearNoise = array2table(awgn(table2array(fear(:, 1:20)),10,0,S));
    fearNoise = [fearNoise, fear(:,21:22)];
    fearNoise.Properties.VariableNames = fear.Properties.VariableNames;
    reset(S);
else
    fearNoise=[];
end
fear = [fear;fearNoise];
%Tabella con Augmentation e Downsampling senza normalizzazione
TAugmented = [happiness; boredom; sadness; anxiety; anger; disgust; fear];

data=TAugmented(:,1:20);
emotion=TAugmented(:,21:22);

results=trainedModel.predictFcn(data);

plotconfusion(categorical(table2array(emotion(:,1))),results);

confronto=[emotion(:,1),array2table(results(:,1))];
