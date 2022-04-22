format longG

sizeData=8064; %numero di misurazioni per canale e traccia
dataRows=40; %numero di traccia
people=32; %numero di partecipanti
windowsNumber = fix(sizeData/(128*10))-1;
dataTrack=zeros(40,people*windowsNumber,20); %tabella in cui metterò nell'ordine traccia,persona,features
labelsTrack=zeros(40,people*windowsNumber,2); %tabella in cui metterò nell'ordine traccia,persona,(valence,arousal)

for n = 1:people  %per ogni partecipante
    disp(strcat("persona: ",num2str(n)));
    matrix=[];
    Responses=[];
    if n<10
        fileName=strcat("data_preprocessed_matlab\s0",num2str(n),".mat");
    else
        fileName=strcat("data_preprocessed_matlab\s",num2str(n),".mat");
    end
    load(fileName); %importa le tabelle data e labels
    
    for i = 1:dataRows  %per ogni traccia
        disp(strcat("traccia: ",num2str(i)));
        
        EEGRAW=zeros(4,sizeData);  %[RAW_TP9,RAW_AF7,RAW_AF8,RAW_TP10]x[8064 campioni]
        EMG=zeros(sizeData,1);
        GSR=zeros(sizeData,1);
        PPG=zeros(sizeData,1);

        %--ESTRAGGO I DATI DEI CANALI CHE MI INTERESSANO---
        for k=1:sizeData
            EEGRAW(1,k)=(data(i,8,k)+data(i,12,k))/2;                              %RAW_TP9
            EEGRAW(2,k)=(data(i,1,k)+data(i,2,k)+data(i,3,k)+data(i,4,k))/4;       %RAW_AF7
            EEGRAW(3,k)=(data(i,17,k)+data(i,18,k)+data(i,20,k)+data(i,21,k))/4;   %RAW_AF8
            EEGRAW(4,k)=(data(i,26,k)+data(i,30,k))/2;                             %RAW_TP10
            EMG(k)=data(i,36,k);
            GSR(k)=data(i,37,k);
            PPG(k)=data(i,39,k);
        end

        %---------CALCOLO FEATURES--------------------------
        for k=0:windowsNumber
            row=[];
            GSRk=GSR((k*windowsNumber*128)+1:(k+1)*windowsNumber*128);%(1+i*windowSize*256:(i+1)*windowSize*256,j)
            PPGk=PPG((k*windowsNumber*128)+1:(k+1)*windowsNumber*128);
            EEGRAWk=EEGRAW(:,(k*windowsNumber*128)+1:(k+1)*windowsNumber*128);
            EMGk=EMG((k*windowsNumber*128)+1:(k+1)*windowsNumber*128);
            
            %calcolo features GSR
            disp(strcat("estrazione feature GSR ",num2str(k)));
            GSRFeatures = [median(GSRk), max(GSRk), rms(GSRk), mean(GSRk, 'omitnan')];
            row=[row,GSRFeatures];
            %calcolo features PPG
            disp(strcat("estrazione feature PPG ",num2str(k)));
            PPGFeatures = [median(PPGk), min(PPGk), max(PPGk)];
            row=[row,PPGFeatures];

            disp(strcat("estrazione feature EEG ",num2str(k)));
            for j = 1:4 %per ogni canale EEG

              % Calcolo features EEG
              if(j==1) %EEG1
                fftEEGRAW = fft(EEGRAWk(j,:));
                betaWave = abs(ifft(bandpass(fftEEGRAW, [13 30], 128)));
                EEGFeatures = [median(EEGRAWk(j,:)), obw(EEGRAWk(j,:)), median(betaWave)];
              elseif(j==2) %EEG2
                fftEEGRAW = fft(EEGRAWk(j,:));
                alphaWave = abs(ifft(bandpass(fftEEGRAW, [7.5 13], 128)));
                EEGFeatures = mean(alphaWave);
              elseif(j==3) %EEG3
                EEGFeatures = median(EEGRAWk(j,:));
              elseif(j==4) %EEG4
                fftEEGRAW = fft(EEGRAWk(j,:));
                thetaWave = abs(ifft(bandpass(fftEEGRAW, [4 7.5], 128)));
                EEGFeatures = [median(EEGRAWk(j,:)), min(thetaWave)];
              end

              row=[row,EEGFeatures];

            end
            %calcolo features EMG
            disp(strcat("estrazione feature EMG ",num2str(k)));
            EMGFeatures = [median(EMGk), max(EMGk), min(EMGk), trimmean(EMGk, 10), mean(EMGk, 'omitnan'), rms(EMGk)];
            row=[row,EMGFeatures]; 

            matrix=[matrix;row];
            dataTrack(i,(n-1)*windowsNumber+k+1,:)=row;
            labelsTrack(i,(n-1)*windowsNumber+k+1,1)=labels(i,1);
            labelsTrack(i,(n-1)*windowsNumber+k+1,2)=labels(i,2);
        end

    end

    Responses=[labels(:,1),labels(:,2)];

    %-----------------NOMINO LE COLONNE------------------------------------
    dataTable = array2table(matrix);
    ResponsesTable = array2table(Responses);
    variableNames={'GSRmedian','GSRmax','GSRrms','GSRmean','PPGmedian','PPGmin','PPGmax','EEG1median','EEG1obw','EEG1betaMedian','EEG2alphaMean','EEG3median','EEG4median','EEG4thetaMin','EMGmedian','EMGmax','EMGmin','EMGtrimmean','EMGmean','EMGrms'}; 
    dataTable.Properties.VariableNames = variableNames;
    ResponsesTable.Properties.VariableNames = {'Valence', 'Arousal'};

    
end

%----------------SALVO LE TABELLE RAGGRUPPATE PER TRACCIA-------------

for k=1:dataRows
  matrix=zeros(people*windowsNumber,20);
  matrix(:,:)=dataTrack(k,1:160,:);
  dataTable= array2table(matrix);
  variableNames={'GSRmedian','GSRmax','GSRrms','GSRmean','PPGmedian','PPGmin','PPGmax','EEG1median','EEG1obw','EEG1betaMedian','EEG2alphaMean','EEG3median','EEG4median','EEG4thetaMin','EMGmedian','EMGmax','EMGmin','EMGtrimmean','EMGmean','EMGrms'}; 
  dataTable.Properties.VariableNames = variableNames;
  writetable(dataTable,strcat("results\Features_t",num2str(k),".csv")); 
  
  matrix=zeros(people*windowsNumber,2);
  matrix(:,:)=labelsTrack(k,1:160,:);
  ResponsesTable = array2table(matrix);
  ResponsesTable.Properties.VariableNames = {'Valence', 'Arousal'};
  writetable(ResponsesTable,strcat("results\Responses_t",num2str(k),".csv"));
end
    




