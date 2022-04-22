format longG

nTracce=40;

for i=1:nTracce
    traccia=readtable(strcat("features/Features_t",num2str(i),".csv"));
    labels=readtable(strcat("features/Responses_t",num2str(i),".csv"));

    Emotion=trainedModel.predictFcn(traccia);
    Valance=labels(:,1);
    Arousal=labels(:,2);
    Emotion=array2table(Emotion);
    Results=[Valance,Arousal,Emotion];

    Results.Properties.VariableNames = {'Valence', 'Arousal','Emotion'};
    writetable(Results,strcat("results/Results_t",num2str(i),".csv"));
end