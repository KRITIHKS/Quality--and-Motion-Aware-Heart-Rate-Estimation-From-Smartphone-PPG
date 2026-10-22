%% Q-MADSP-HR: subject-clustered comparison and reliability ablation
% Uses the actual column names in record_level_results.csv.
% Run in a folder containing record_level_results.csv.
% This script analyzes saved predictions; it does NOT regenerate predictions.
clear; clc;
ledgerFile = 'record_level_results.csv';
if ~isfile(ledgerFile)
    error('Place record_level_results.csv in the current folder.');
end
T = readtable(ledgerFile,'VariableNamingRule','preserve');
required = {'SubjectID','ReferenceHR','PeakHR','AdaptiveHR','Reliability', ...
    'EstimatorAgreement','Periodicity','SpectralConcentration', ...
    'SpectralEntropy','MotionIndex','Status'};
missing = setdiff(required,T.Properties.VariableNames);
if ~isempty(missing), error('Missing columns: %s',strjoin(missing,', ')); end

T = T(strcmpi(string(T.Status),'SUCCESS') & isfinite(T.ReferenceHR) & ...
    isfinite(T.AdaptiveHR) & ~ismissing(string(T.SubjectID)),:);
T.SubjectID = string(T.SubjectID);
subs = unique(T.SubjectID);
rng(2026,'twister');

% Paired adaptive-vs-peak comparison; only matched finite estimates.
D = T(isfinite(T.PeakHR),:);
D.Delta = abs(D.AdaptiveHR-D.ReferenceHR)-abs(D.PeakHR-D.ReferenceHR);
obsDelta = mean(D.Delta);
B = 5000; bootDelta = nan(B,1);
dSubs = unique(D.SubjectID);
for b=1:B
    picked = dSubs(randi(numel(dSubs),numel(dSubs),1));
    vals = [];
    for j=1:numel(picked)
        vals = [vals; D.Delta(D.SubjectID==picked(j))]; %#ok<AGROW>
    end
    bootDelta(b)=mean(vals);
end
ciDelta = prctile(bootDelta,[2.5 97.5]);
% Compute subject-specific MAE correctly from absolute errors.
D.AdaptiveAbsError = abs(D.AdaptiveHR-D.ReferenceHR);
D.PeakAbsError = abs(D.PeakHR-D.ReferenceHR);
subj = groupsummary(D,'SubjectID','mean',{'AdaptiveAbsError','PeakAbsError'});
[pWil,~,statsWil] = signrank(subj.mean_AdaptiveAbsError,subj.mean_PeakAbsError);
paired = table(height(D),numel(dSubs),mean(D.AdaptiveAbsError),mean(D.PeakAbsError), ...
    obsDelta,ciDelta(1),ciDelta(2),pWil, ...
    'VariableNames',{'Npaired','Nsubjects','AdaptiveMAE_bpm','PeakMAE_bpm', ...
    'DeltaMAE_AdaptiveMinusPeak','CI95Lower','CI95Upper','SubjectLevelWilcoxonP'});
writetable(paired,'Q_MADSP_HR_paired_subject_cluster_results.csv');

% Subject-cluster bootstrap AUC for classifying |adaptive error| <= 10 bpm.
T.AbsError = abs(T.AdaptiveHR-T.ReferenceHR);
T.Good10 = T.AbsError<=10;
scoreNames = {'Reliability','EstimatorAgreement','Periodicity', ...
    'SpectralConcentration','SpectralEntropy','MotionIndex'};
scoreLabels = {'Composite reliability','Estimator agreement','Periodicity', ...
    'Spectral concentration','Inverse spectral entropy','Inverse motion'};
orient = [1 1 1 1 -1 -1]; % higher oriented score means more likely |error|<=10
aucFcn = @(y,s) localAUC(y,s);
pointAUC = nan(numel(scoreNames),1);
for j=1:numel(scoreNames)
    pointAUC(j)=aucFcn(T.Good10,orient(j)*T.(scoreNames{j}));
end
Bauc=3000; bootAUC=nan(Bauc,numel(scoreNames));
for b=1:Bauc
    picked = subs(randi(numel(subs),numel(subs),1));
    ix = false(height(T),1);
    Tb = T([],:);
    for j=1:numel(picked)
        Tb=[Tb; T(T.SubjectID==picked(j),:)]; %#ok<AGROW>
    end
    for j=1:numel(scoreNames)
        bootAUC(b,j)=aucFcn(Tb.Good10,orient(j)*Tb.(scoreNames{j}));
    end
end
lo=nan(numel(scoreNames),1); hi=lo;
for j=1:numel(scoreNames)
    lo(j)=prctile(bootAUC(:,j),2.5);
    hi(j)=prctile(bootAUC(:,j),97.5);
end
ablation=table(string(scoreLabels(:)),pointAUC,lo,hi, ...
    'VariableNames',{'Score','AUC_for_AbsoluteErrorLE10','ClusterBootstrap95Lower','ClusterBootstrap95Upper'});
writetable(ablation,'Q_MADSP_HR_reliability_ablation_results.csv');
deltaAUC=bootAUC(:,1)-bootAUC(:,2);
aucDiff=table(pointAUC(1)-pointAUC(2),prctile(deltaAUC,2.5),prctile(deltaAUC,97.5), ...
    'VariableNames',{'AUC_Difference_ReliabilityMinusAgreement','ClusterBootstrap95Lower','ClusterBootstrap95Upper'});
writetable(aucDiff,'Q_MADSP_HR_reliability_AUC_difference.csv');

% Descriptive Spearman correlations; do not interpret recording-level p-values
% as population-level inference because recordings are clustered by subject.
rho=nan(numel(scoreNames),1); pdesc=rho;
for j=1:numel(scoreNames)
    [rho(j),pdesc(j)] = corr(orient(j)*T.(scoreNames{j}),T.AbsError,'Type','Spearman','Rows','complete');
end
corrTable=table(string(scoreLabels(:)),rho,pdesc, ...
    'VariableNames',{'Score','SpearmanRhoWithAbsoluteError','RecordingLevelP_Descriptive'});
writetable(corrTable,'Q_MADSP_HR_score_error_correlations.csv');

fprintf('Paired adaptive-minus-peak MAE delta: %.3f bpm (cluster bootstrap 95%% CI %.3f to %.3f); subject-level Wilcoxon p=%.3f\\n', ...
    obsDelta,ciDelta(1),ciDelta(2),pWil);
fprintf('Reliability AUC %.3f; agreement AUC %.3f; AUC difference %.3f (95%% cluster CI %.3f to %.3f)\\n', ...
    pointAUC(1),pointAUC(2),pointAUC(1)-pointAUC(2),prctile(deltaAUC,2.5),prctile(deltaAUC,97.5));

function a = localAUC(y,s)
y = logical(y(:)); s=s(:);
ok=isfinite(s) & ~isnan(double(y));
y=y(ok); s=s(ok);
if numel(unique(y))<2, a=NaN; return; end
[~,~,~,a]=perfcurve(y,s,true);
end
