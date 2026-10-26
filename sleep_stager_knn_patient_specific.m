%% KNN PATIENT-SPECIFIC SLEEP STAGING (with a lot of data leakage!)
% Previously calculated inputs:
% all_feat_eeg  : N x 12
% wind_label_eeg: 1 x N or N x 1
%
% Labels: 0=W, 1=N1, 2=N2, 3=N3, 4=REM, 9=Extra

clearvars -except all_feat_eeg wind_label_eeg
clc; close all;
rng(42);

%% 1. Data preparation

X = double(all_feat_eeg);
labels = double(wind_label_eeg(:));

assert(size(X,1) == numel(labels), ...
    'Number of features and labels does not match.');

% Remove extra labels, NaN values, and non-finite features
valid = ismember(labels,[0 1 2 3 4]) & ...
        all(isfinite(X),2);

X = X(valid,:);
labels = labels(valid);

% Fixed class order
classValues = [0 1 2 3 4];
classNames = {'W','N1','N2','N3','REM'};

Y = categorical(labels, classValues, classNames);

fprintf('Valid windows: %d\n',size(X,1));
disp(countcats(Y));

%% 2. Chronological split
% Patient-specific: all windows belong to the same patient.
% The windows must be in chronological order.
% With 50% overlap, each window shares samples with adjacent windows.

N = size(X,1);
assert(N >= 50,'Insufficient number of windows.');

i1 = floor(0.60*N);
i2 = floor(0.80*N);

% Exclude 2 windows around the boundaries
gap = 2;

idxTrain = (1:i1-gap)';
idxVal   = (i1+gap+1:i2-gap)';
idxTest  = (i2+gap+1:N)';

assert(~isempty(idxTrain) && ~isempty(idxVal) && ...
       ~isempty(idxTest),'Empty partitions.');

XTrain = X(idxTrain,:);
YTrain = Y(idxTrain);

XVal = X(idxVal,:);
YVal = Y(idxVal);

XTest = X(idxTest,:);
YTest = Y(idxTest);

fprintf('Training:   %d windows\n',numel(idxTrain));
fprintf('Validation: %d windows\n',numel(idxVal));
fprintf('Test:       %d windows\n',numel(idxTest));

% Check that the training set contains all expected classes
present = categories(removecats(YTrain));
disp('Classes in the training set:');
disp(present);

%% 3. Select K using the validation set

Klist = [1 3 5 7 9 11 15 21 31];
Klist = Klist(Klist < size(XTrain,1));

macroF1Val = nan(size(Klist));
accuracyVal = nan(size(Klist));

for q = 1:numel(Klist)

    mdlTmp = fitcknn(XTrain,YTrain, ...
        'NumNeighbors',Klist(q), ...
        'Distance','euclidean', ...
        'DistanceWeight','inverse', ...
        'Standardize',true, ...
        'ClassNames',categorical(classNames,classNames));

    predVal = predict(mdlTmp,XVal);

    [macroF1Val(q), accuracyVal(q)] = ...
        getMetrics(YVal,predVal,classNames);

    fprintf('K=%2d | Accuracy=%.3f | Macro-F1=%.3f\n', ...
        Klist(q),accuracyVal(q),macroF1Val(q));
end

% Select K based on macro-F1, not accuracy alone
[~,bestIdx] = max(macroF1Val);
bestK = Klist(bestIdx);

fprintf('\nOptimal K: %d\n',bestK);

%% 4. Validation plot

figure('Name','KNN validation','Color','w');

plot(Klist,accuracyVal,'-o','LineWidth',1.8);
hold on;
plot(Klist,macroF1Val,'-s','LineWidth',1.8);
xline(bestK,'--k',sprintf('K=%d',bestK));

xlabel('Number of neighbors K');
ylabel('Performance');
title('K selection using the validation set');
legend('Accuracy','Macro-F1','Location','best');
grid on;

%% 5. Final model
% Refit on training + validation after selecting K.
% The test set remains completely separate.

XFit = [XTrain; XVal];
YFit = [YTrain; YVal];

knnModel = fitcknn(XFit,YFit, ...
    'NumNeighbors',bestK, ...
    'Distance','euclidean', ...
    'DistanceWeight','inverse', ...
    'Standardize',true, ...
    'ClassNames',categorical(classNames,classNames));

YPred = predict(knnModel,XTest);

%% 6. Test set metrics

[macroF1Test,accuracyTest,precision,recall,f1,C] = ...
    getMetrics(YTest,YPred,classNames);

fprintf('\n===== TEST RESULTS =====\n');
fprintf('Accuracy : %.4f (%.2f%%)\n', ...
    accuracyTest,100*accuracyTest);
fprintf('Macro-F1 : %.4f\n',macroF1Test);

Tmetrics = table(string(classNames(:)),precision,recall,f1, ...
    'VariableNames',{'Stage','Precision','Recall','F1'});
disp(Tmetrics);

%% 7. Confusion matrix (MATLAB compatibility)

classNames = {'W','N1','N2','N3','REM'};

% Fixed class order
classOrder = categorical(1:5, 1:5, classNames);

YTestCat = categorical(double(YTest), 1:5, classNames);
YPredCat = categorical(double(YPred), 1:5, classNames);

% Explicitly calculate the confusion matrix
C = confusionmat(YTestCat, YPredCat, ...
    'Order', classOrder);

% Display the matrix without passing the Order property
figure('Name','Confusion matrix','Color','w');

cm = confusionchart(C, classNames);

cm.Title = sprintf('KNN Patient-Specific | K=%d', bestK);
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';

%% 8. Per-class performance

figure('Name','Per-class metrics','Color','w');

bar([precision recall f1]);
set(gca,'XTick',1:numel(classNames), ...
    'XTickLabel',classNames);

ylabel('Score');
ylim([0 1]);
title('Metrics by sleep stage - Test set');
legend('Precision','Recall','F1','Location','best');
grid on;

%% 9. Temporal sequence: ground-truth vs. predicted labels

% Convert sleep stages to numerical codes 1-5
yTrueNum = double(YTest);
yPredNum = double(YPred);

figure('Name','Hypnogram comparison','Color','w');

subplot(2,1,1);
stairs(idxTest,yTrueNum,'LineWidth',1.2);
yticks(1:5);
yticklabels(classNames);
ylim([0.5 5.5]);
title('Ground-truth sleep stages - Test set');
ylabel('Sleep stage');
grid on;

subplot(2,1,2);
stairs(idxTest,yPredNum,'LineWidth',1.2);
yticks(1:5);
yticklabels(classNames);
ylim([0.5 5.5]);
title('Predicted sleep stages - KNN');
xlabel('Window index');
ylabel('Sleep stage');
grid on;

% figure('Name','Hypnogram comparison','Color','w');
%
% subplot(2,1,1);
% stairs(idxTest,YTest,'LineWidth',1.2);
% yticks(1:5);
% yticklabels(classNames);
% ylim([0.5 5.5]);
% title('Ground-truth sleep stages - Test set');
% ylabel('Sleep stage');
% grid on;
%
% subplot(2,1,2);
% stairs(idxTest,YPred,'LineWidth',1.2);
% yticks(1:5);
% yticklabels(classNames);
% ylim([0.5 5.5]);
% title('Predicted sleep stages - KNN');
% xlabel('Window index');
% ylabel('Sleep stage');
% grid on;

%% Local function: multiclass metrics
function [macroF1,accuracy,precision,recall,f1,C] = ...
    getMetrics(Ytrue,Ypred,classNames)

    order = categorical(classNames,classNames);

    Ytrue = categorical(string(Ytrue),classNames);
    Ypred = categorical(string(Ypred),classNames);

    C = confusionmat(Ytrue,Ypred,'Order',order);

    TP = diag(C);
    precision = TP ./ max(sum(C,1)',1);
    recall = TP ./ max(sum(C,2),1);

    f1 = 2*precision.*recall ./ ...
        max(precision+recall,eps);

    accuracy = sum(TP)/max(sum(C(:)),1);
    macroF1 = mean(f1);
end
