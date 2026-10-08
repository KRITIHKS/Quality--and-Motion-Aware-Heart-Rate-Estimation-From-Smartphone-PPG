Y
%% ================================================================
%  BUT PPG v2.0.0
%  QUALITY-DRIVEN ADAPTIVE DSP FOR SMARTPHONE HEART-RATE ESTIMATION
%
%  JOURNAL-ORIENTED EXPERIMENTAL PIPELINE
%
%  Major features:
%    1. Robust local WFDB loading
%    2. Metadata verification
%    3. Filesystem integrity check
%    4. Subject-independent 5-fold partition
%    5. Five candidate HR estimators
%         - Peak
%         - FFT
%         - Welch PSD
%         - ACF
%         - Harmonic
%    6. Reference-free signal-quality features
%    7. Motion contamination analysis
%    8. Adaptive reliability-aware estimator selection
%    9. Selective/risk-coverage analysis
%   10. Matched-coverage comparison
%   11. Quality / motion / HR stratification
%   12. Subject-level analysis
%   13. Subject-level bootstrap confidence intervals
%   14. Paired statistical testing
%   15. Bland-Altman analysis
%   16. Publication figures and CSV/MAT outputs
%
%  IMPORTANT:
%    Results are computed from the actual signals.
%    No target MAE/accuracy is imposed.
%
% ================================================================

clear;
clc;
close all;

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' BUT PPG JOURNAL V2 - FIXED LOCAL WFDB PIPELINE\n');
fprintf(' Quality-Driven Adaptive DSP Heart-Rate Estimation\n');
fprintf('===============================================================\n\n');

%% ================================================================
% 1. USER CONFIGURATION
% ================================================================

DATA_ROOT = ...
    'C:\Users\ASUS\Downloads\brno-university-of-technology-smartphone-ppg-database-but-ppg-2.0.0\brno-university-of-technology-smartphone-ppg-database-but-ppg-2.0.0';

WFDB_PATH = 'C:\WFDB\mcode';

RESULT_ROOT = fullfile(DATA_ROOT,'BUT_PPG_JOURNAL_V2_FIXED_RESULTS');

FIG_ROOT = fullfile(RESULT_ROOT,'figures');

if ~exist(RESULT_ROOT,'dir')
    mkdir(RESULT_ROOT);
end

if ~exist(FIG_ROOT,'dir')
    mkdir(FIG_ROOT);
end

%% Signal configuration

FS_TARGET = 30;

LOW_HZ  = 0.7;
HIGH_HZ = 3.5;

MIN_HR = 30;
MAX_HR = 180;

MIN_SECONDS = 8;

NFFT = 4096;

%% Evaluation configuration

N_FOLDS = 5;

RNG_SEED = 2026;

TARGET_COVERAGE = 0.40;

N_BOOT = 2000;

CHECKPOINT_EVERY = 100;

%% ================================================================
% 2. MATLAB / WFDB SETUP
% ================================================================

fprintf('Checking WFDB installation...\n');

if exist(WFDB_PATH,'dir')
    addpath(genpath(WFDB_PATH));
    fprintf('WFDB path added:\n%s\n',WFDB_PATH);
else
    warning('WFDB directory not found: %s',WFDB_PATH);
end

fprintf('\nDataset root:\n%s\n',DATA_ROOT);

if ~exist(DATA_ROOT,'dir')
    error('DATA_ROOT does not exist.');
end

%% ================================================================
% 3. FIND METADATA FILES
% ================================================================

fprintf('\nLoading metadata...\n');

qualityCandidates = { ...
    'quality-hr-ann.csv', ...
    'quality-hr-ann(2).csv'};

subjectCandidates = { ...
    'subject-info.csv', ...
    'subject-info(2).csv'};

qualityFile = '';

for k = 1:numel(qualityCandidates)
    candidate = fullfile(DATA_ROOT,qualityCandidates{k});
    if exist(candidate,'file')
        qualityFile = candidate;
        break;
    end
end

if isempty(qualityFile)
    error('Could not find quality-hr-ann.csv in DATA_ROOT.');
end

subjectFile = '';

for k = 1:numel(subjectCandidates)
    candidate = fullfile(DATA_ROOT,subjectCandidates{k});
    if exist(candidate,'file')
        subjectFile = candidate;
        break;
    end
end

if isempty(subjectFile)
    warning('subject-info.csv not found. Subject metadata will be limited.');
end

fprintf('Quality metadata:\n%s\n',qualityFile);

if ~isempty(subjectFile)
    fprintf('Subject metadata:\n%s\n',subjectFile);
end

%% ================================================================
% 4. READ QUALITY METADATA
% ================================================================

Q = readtable(qualityFile,'VariableNamingRule','preserve');

fprintf('\nQuality metadata columns:\n');
disp(Q.Properties.VariableNames');

%% Locate columns robustly

idCol = findColumn(Q,{'ID','id','Signal ID','signal_id'});

qualityCol = findColumn(Q,{'Quality','quality'});

hrCol = findColumn(Q,{'HR','hr','Reference HR','reference_hr'});

if isempty(idCol)
    error('ID column not found in quality metadata.');
end

if isempty(qualityCol)
    error('Quality column not found in quality metadata.');
end

if isempty(hrCol)
    error('HR column not found in quality metadata.');
end

ID = double(Q.(idCol));

Quality = double(Q.(qualityCol));

ReferenceHR = double(Q.(hrCol));

%% Remove invalid metadata

validMeta = ...
    isfinite(ID) & ...
    isfinite(Quality) & ...
    isfinite(ReferenceHR);

ID = ID(validMeta);
Quality = Quality(validMeta);
ReferenceHR = ReferenceHR(validMeta);

fprintf('\nValid metadata records: %d\n',numel(ID));

%% ================================================================
% 5. READ SUBJECT METADATA
% ================================================================

SubjectTable = table();

if ~isempty(subjectFile)

    S = readtable(subjectFile,'VariableNamingRule','preserve');

    sidCol = findColumn(S,{'ID','id','Signal ID','signal_id'});

    if ~isempty(sidCol)

        sid = double(S.(sidCol));

        [~,ia,ib] = intersect(ID,sid);

        SubjectTable = S(ib,:);

        % Reorder according to quality table IDs
        temp = nan(numel(ID),1);

        temp(ia) = 1:numel(ia);

        [~,order] = sort(temp);

        if all(isfinite(temp))
            SubjectTable = SubjectTable(order,:);
        end

    end
end

%% ================================================================
% 6. DERIVE SUBJECT IDS
% ================================================================

% BUT PPG record IDs are six-digit IDs.
% First three digits identify the subject.

SubjectID = floor(ID/1000);

fprintf('Unique subjects: %d\n',numel(unique(SubjectID)));

fprintf('Reference HR range: %.1f - %.1f BPM\n', ...
    min(ReferenceHR),max(ReferenceHR));

fprintf('Quality = 1 records: %d\n',sum(Quality==1));
fprintf('Quality = 0 records: %d\n',sum(Quality==0));

%% ================================================================
% 7. VERIFY DATASET FILESYSTEM
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' FILESYSTEM INTEGRITY CHECK\n');
fprintf('===============================================================\n');

N = numel(ID);

PPGExists = false(N,1);
ACCExists = false(N,1);

for i = 1:N

    rid = sprintf('%06d',round(ID(i)));

    recordDir = fullfile(DATA_ROOT,rid);

    ppgBase = fullfile(recordDir,[rid '_PPG']);
    accBase = fullfile(recordDir,[rid '_ACC']);

    PPGExists(i) = ...
        exist([ppgBase '.hea'],'file') && ...
        exist([ppgBase '.dat'],'file');

    ACCExists(i) = ...
        exist([accBase '.hea'],'file') && ...
        exist([accBase '.dat'],'file');

end

fprintf('Total metadata records : %d\n',N);
fprintf('PPG files found        : %d/%d (%.2f%%)\n', ...
    sum(PPGExists),N,100*mean(PPGExists));

fprintf('ACC files found        : %d/%d (%.2f%%)\n', ...
    sum(ACCExists),N,100*mean(ACCExists));

if sum(PPGExists)==0
    error('No PPG records found. DATA_ROOT/path is incorrect.');
end

%% ================================================================
% 8. PRE-FLIGHT WFDB TEST
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' WFDB LOCAL LOADER PRE-FLIGHT TEST\n');
fprintf('===============================================================\n');

firstPPG = find(PPGExists,1,'first');

testRID = sprintf('%06d',round(ID(firstPPG)));

testPPGBase = fullfile( ...
    DATA_ROOT, ...
    testRID, ...
    [testRID '_PPG']);

fprintf('Testing record:\n%s\n',testPPGBase);

[testSignal,testFs,testOK,testMessage] = ...
    loadLocalWFDBRecord(testPPGBase);

if ~testOK

    fprintf('\nWFDB TEST FAILED.\n');
    fprintf('%s\n',testMessage);

    error(['Local WFDB loading failed. ', ...
           'The script stopped before processing the dataset.']);
end

fprintf('WFDB test successful.\n');
fprintf('Samples : %d\n',size(testSignal,1));
fprintf('Fs      : %.3f Hz\n',testFs);

clear testSignal testFs;

%% ================================================================
% 9. SUBJECT-INDEPENDENT 5-FOLD ASSIGNMENT
% ================================================================

fprintf('\nCreating subject-independent folds...\n');

rng(RNG_SEED);

uniqueSubjects = unique(SubjectID);

uniqueSubjects = uniqueSubjects(randperm(numel(uniqueSubjects)));

FoldID = zeros(N,1);

for s = 1:numel(uniqueSubjects)

    subj = uniqueSubjects(s);

    fold = mod(s-1,N_FOLDS)+1;

    FoldID(SubjectID==subj) = fold;

end

fprintf('Number of folds: %d\n',N_FOLDS);

for f = 1:N_FOLDS

    fprintf(' Fold %d: %d records, %d subjects\n', ...
        f, ...
        sum(FoldID==f), ...
        numel(unique(SubjectID(FoldID==f))));

end

%% ================================================================
% 10. PREALLOCATE RESULT TABLE
% ================================================================

MethodNames = { ...
    'Peak', ...
    'FFT', ...
    'PSD', ...
    'ACF', ...
    'Harmonic', ...
    'Adaptive'};

M = numel(MethodNames);

R = table();

R.ID = ID;
R.SubjectID = SubjectID;
R.Fold = FoldID;

R.Quality = Quality;
R.ReferenceHR = ReferenceHR;

R.PPGExists = PPGExists;
R.ACCExists = ACCExists;

R.PeakHR = nan(N,1);
R.FFTHR = nan(N,1);
R.PSDHR = nan(N,1);
R.ACFHR = nan(N,1);
R.HarmonicHR = nan(N,1);
R.AdaptiveHR = nan(N,1);

%% Quality features

R.Periodicity = nan(N,1);
R.SpectralConcentration = nan(N,1);
R.SpectralEntropy = nan(N,1);
R.EstimatorAgreement = nan(N,1);
R.MotionIndex = nan(N,1);

R.Reliability = nan(N,1);

R.SelectedMethod = strings(N,1);

R.Status = strings(N,1);
R.ErrorMessage = strings(N,1);

%% ================================================================
% 11. PROCESS ALL RECORDS
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' PROCESSING RECORDS\n');
fprintf('===============================================================\n');

checkpointFile = fullfile(RESULT_ROOT,'checkpoint.mat');

startIndex = 1;

%% Resume if checkpoint exists

if exist(checkpointFile,'file')

    fprintf('\nExisting checkpoint detected.\n');

    answer = input('Resume from checkpoint? Y/N [Y]: ','s');

    if isempty(answer)
        answer = 'Y';
    end

    if upper(answer(1))=='Y'

        C = load(checkpointFile);

        if isfield(C,'R') && height(C.R)==N

            R = C.R;

            done = R.Status ~= "";

            if any(done)

                lastDone = find(done,1,'last');

                startIndex = lastDone + 1;

                fprintf('Resuming from record %d.\n',startIndex);

            end

        end
    end
end

tic;

for i = startIndex:N

    rid = sprintf('%06d',round(ID(i)));

    if ~PPGExists(i)

        R.Status(i) = "SKIPPED_NO_PPG";
        R.ErrorMessage(i) = "PPG files not found.";

        continue;
    end

    try

        [out,ok,msg] = processOneRecord( ...
            DATA_ROOT, ...
            rid, ...
            FS_TARGET, ...
            LOW_HZ, ...
            HIGH_HZ, ...
            MIN_HR, ...
            MAX_HR, ...
            MIN_SECONDS, ...
            NFFT);

        if ok

            R.PeakHR(i) = out.PeakHR;
            R.FFTHR(i) = out.FFTHR;
            R.PSDHR(i) = out.PSDHR;
            R.ACFHR(i) = out.ACFHR;
            R.HarmonicHR(i) = out.HarmonicHR;
            R.AdaptiveHR(i) = out.AdaptiveHR;

            R.Periodicity(i) = out.Periodicity;
            R.SpectralConcentration(i) = out.SpectralConcentration;
            R.SpectralEntropy(i) = out.SpectralEntropy;
            R.EstimatorAgreement(i) = out.EstimatorAgreement;
            R.MotionIndex(i) = out.MotionIndex;
            R.Reliability(i) = out.Reliability;

            R.SelectedMethod(i) = out.SelectedMethod;

            R.Status(i) = "SUCCESS";
            R.ErrorMessage(i) = "";

        else

            R.Status(i) = "FAILED";
            R.ErrorMessage(i) = string(msg);

            fprintf('Record %s failed: %s\n',rid,msg);

        end

    catch ME

        R.Status(i) = "FAILED";
        R.ErrorMessage(i) = string(ME.message);

        fprintf('Record %s failed: %s\n', ...
            rid,ME.message);

    end

    %% Progress

    if mod(i,CHECKPOINT_EVERY)==0 || i==N

        elapsed = toc;

        processed = i-startIndex+1;

        rate = processed/max(elapsed,eps);

        remaining = N-i;

        etaMin = remaining/max(rate,eps)/60;

        fprintf('\nProgress: %d/%d (%.1f%%)\n', ...
            i,N,100*i/N);

        fprintf('Elapsed: %.1f min | ETA: %.1f min\n', ...
            elapsed/60,etaMin);

        save(checkpointFile,'R','-v7.3');

        fprintf('Checkpoint saved.\n');

    end

end

%% ================================================================
% 12. BASIC PROCESSING SUMMARY
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' PROCESSING SUMMARY\n');
fprintf('===============================================================\n');

success = R.Status=="SUCCESS";

fprintf('Total records      : %d\n',height(R));
fprintf('Candidate-ready    : %d (%.2f%%)\n', ...
    sum(PPGExists),100*mean(PPGExists));

fprintf('Successful         : %d (%.2f%%)\n', ...
    sum(success),100*mean(success));

fprintf('Failed/skipped     : %d\n',sum(~success));

%% Save raw results

writetable(R,fullfile(RESULT_ROOT,'record_level_results.csv'));

save(fullfile(RESULT_ROOT,'record_level_results.mat'), ...
    'R','-v7.3');

%% ================================================================
% 13. VALIDITY MASKS
% ================================================================

reference = R.ReferenceHR;

validAdaptive = ...
    success & ...
    isfinite(R.AdaptiveHR) & ...
    isfinite(reference);

validPeak = ...
    success & ...
    isfinite(R.PeakHR) & ...
    isfinite(reference);

validFFT = ...
    success & ...
    isfinite(R.FFTHR) & ...
    isfinite(reference);

validPSD = ...
    success & ...
    isfinite(R.PSDHR) & ...
    isfinite(reference);

validACF = ...
    success & ...
    isfinite(R.ACFHR) & ...
    isfinite(reference);

validHarmonic = ...
    success & ...
    isfinite(R.HarmonicHR) & ...
    isfinite(reference);

%% ================================================================
% 14. OVERALL PERFORMANCE
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' OVERALL PERFORMANCE\n');
fprintf('===============================================================\n');

metricTable = table();

for m = 1:M

    method = MethodNames{m};

    switch method

        case 'Peak'
            pred = R.PeakHR;
            mask = validPeak;

        case 'FFT'
            pred = R.FFTHR;
            mask = validFFT;

        case 'PSD'
            pred = R.PSDHR;
            mask = validPSD;

        case 'ACF'
            pred = R.ACFHR;
            mask = validACF;

        case 'Harmonic'
            pred = R.HarmonicHR;
            mask = validHarmonic;

        case 'Adaptive'
            pred = R.AdaptiveHR;
            mask = validAdaptive;

    end

    met = calculateMetrics( ...
        reference(mask), ...
        pred(mask));

    row = table( ...
        string(method), ...
        sum(mask), ...
        met.MAE, ...
        met.RMSE, ...
        met.MAPE, ...
        met.Bias, ...
        met.AER5, ...
        met.AER10, ...
        met.PearsonR, ...
        met.R2, ...
        'VariableNames',{ ...
        'Method','N','MAE','RMSE','MAPE','Bias', ...
        'AER5','AER10','PearsonR','R2'});

    metricTable = [metricTable;row];

end

disp(metricTable);

writetable(metricTable, ...
    fullfile(RESULT_ROOT,'overall_metrics.csv'));

%% ================================================================
% 15. 5-FOLD SUBJECT-INDEPENDENT RESULTS
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' SUBJECT-INDEPENDENT 5-FOLD RESULTS\n');
fprintf('===============================================================\n');

FoldResults = table();

for f = 1:N_FOLDS

    maskFold = validAdaptive & R.Fold==f;

    if sum(maskFold)<2
        continue;
    end

    met = calculateMetrics( ...
        reference(maskFold), ...
        R.AdaptiveHR(maskFold));

    row = table( ...
        f, ...
        sum(maskFold), ...
        numel(unique(R.SubjectID(maskFold))), ...
        met.MAE, ...
        met.RMSE, ...
        met.MAPE, ...
        met.Bias, ...
        met.AER5, ...
        met.AER10, ...
        met.PearsonR, ...
        met.R2, ...
        'VariableNames',{ ...
        'Fold','N','Subjects','MAE','RMSE','MAPE','Bias', ...
        'AER5','AER10','PearsonR','R2'});

    FoldResults = [FoldResults;row];

end

disp(FoldResults);

writetable(FoldResults, ...
    fullfile(RESULT_ROOT,'five_fold_results.csv'));

%% ================================================================
% 16. SUBJECT-LEVEL PERFORMANCE
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' SUBJECT-LEVEL PERFORMANCE\n');
fprintf('===============================================================\n');

subjects = unique(R.SubjectID(validAdaptive));

SubjectResults = table();

for s = 1:numel(subjects)

    sid = subjects(s);

    mask = validAdaptive & R.SubjectID==sid;

    if sum(mask)<1
        continue;
    end

    met = calculateMetrics( ...
        reference(mask), ...
        R.AdaptiveHR(mask));

    row = table( ...
        sid, ...
        sum(mask), ...
        met.MAE, ...
        met.RMSE, ...
        met.MAPE, ...
        met.Bias, ...
        met.AER5, ...
        met.AER10, ...
        met.PearsonR, ...
        met.R2, ...
        'VariableNames',{ ...
        'SubjectID','N','MAE','RMSE','MAPE','Bias', ...
        'AER5','AER10','PearsonR','R2'});

    SubjectResults = [SubjectResults;row];

end

disp(SubjectResults);

writetable(SubjectResults, ...
    fullfile(RESULT_ROOT,'subject_level_results.csv'));

%% ================================================================
% 17. BOOTSTRAP CONFIDENCE INTERVALS
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' SUBJECT-LEVEL BOOTSTRAP 95%% CIs\n');
fprintf('===============================================================\n');

BootstrapResults = table();

for m = 1:M

    method = MethodNames{m};

    switch method

        case 'Peak'
            pred = R.PeakHR;
            mask = validPeak;

        case 'FFT'
            pred = R.FFTHR;
            mask = validFFT;

        case 'PSD'
            pred = R.PSDHR;
            mask = validPSD;

        case 'ACF'
            pred = R.ACFHR;
            mask = validACF;

        case 'Harmonic'
            pred = R.HarmonicHR;
            mask = validHarmonic;

        case 'Adaptive'
            pred = R.AdaptiveHR;
            mask = validAdaptive;

    end

    [bootMAE,ciMAE] = ...
        subjectBootstrapMAE( ...
        R.SubjectID(mask), ...
        reference(mask), ...
        pred(mask), ...
        N_BOOT,RNG_SEED);

    row = table( ...
        string(method), ...
        bootMAE, ...
        ciMAE(1), ...
        ciMAE(2), ...
        'VariableNames',{ ...
        'Method','BootstrapMAE','CI95_Low','CI95_High'});

    BootstrapResults = [BootstrapResults;row];

end

disp(BootstrapResults);

writetable(BootstrapResults, ...
    fullfile(RESULT_ROOT,'bootstrap_CI_results.csv'));

%% ================================================================
% 18. RISK-COVERAGE ANALYSIS
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' RISK-COVERAGE ANALYSIS\n');
fprintf('===============================================================\n');

RC = riskCoverageTable(R);

disp(RC);

writetable(RC, ...
    fullfile(RESULT_ROOT,'risk_coverage.csv'));

%% ================================================================
% 19. MATCHED-COVERAGE COMPARISON
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' MATCHED-COVERAGE COMPARISON\n');
fprintf('===============================================================\n');

MC = matchedCoverageTable( ...
    R, ...
    TARGET_COVERAGE);

disp(MC);

writetable(MC, ...
    fullfile(RESULT_ROOT,'matched_coverage.csv'));

%% ================================================================
% 20. QUALITY STRATIFICATION
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' QUALITY STRATIFICATION\n');
fprintf('===============================================================\n');

qualityResults = stratifiedMetrics( ...
    R.Quality, ...
    [0 0.5 1.5], ...
    {'Poor','Good'}, ...
    R.AdaptiveHR, ...
    reference, ...
    validAdaptive);

disp(qualityResults);

writetable(qualityResults, ...
    fullfile(RESULT_ROOT,'quality_stratification.csv'));

%% ================================================================
% 21. HR STRATIFICATION
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' HEART-RATE STRATIFICATION\n');
fprintf('===============================================================\n');

hrResults = stratifiedMetrics( ...
    reference, ...
    [0 60 100 140 250], ...
    {'Low','Normal','High','VeryHigh'}, ...
    R.AdaptiveHR, ...
    reference, ...
    validAdaptive);

disp(hrResults);

writetable(hrResults, ...
    fullfile(RESULT_ROOT,'hr_stratification.csv'));

%% ================================================================
% 22. MOTION STRATIFICATION
% ================================================================

motionValid = validAdaptive & isfinite(R.MotionIndex);

if sum(motionValid)>10

    motionVals = R.MotionIndex(motionValid);

    q1 = quantile(motionVals,1/3);
    q2 = quantile(motionVals,2/3);

    fprintf('\nMotion thresholds: %.4f / %.4f\n',q1,q2);

    motionResults = stratifiedMetrics( ...
        R.MotionIndex, ...
        [-inf q1 q2 inf], ...
        {'LowMotion','MediumMotion','HighMotion'}, ...
        R.AdaptiveHR, ...
        reference, ...
        validAdaptive & motionValid);

    disp(motionResults);

    writetable(motionResults, ...
        fullfile(RESULT_ROOT,'motion_stratification.csv'));

else

    motionResults = table();

end

%% ================================================================
% 23. PAIRED COMPARISON
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' PAIRED METHOD COMPARISON\n');
fprintf('===============================================================\n');

ComparisonResults = pairedComparisonTable(R);

disp(ComparisonResults);

writetable(ComparisonResults, ...
    fullfile(RESULT_ROOT,'paired_comparison.csv'));

%% ================================================================
% 24. FIGURE 1 - METHOD MAE
% ================================================================

figure('Color','w');

bar(categorical(metricTable.Method),metricTable.MAE);

ylabel('MAE (BPM)');
xlabel('Estimator');
title('Heart-Rate Estimation Error by Method');

grid on;

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig1_Method_MAE.png'));

%% ================================================================
% 25. FIGURE 2 - METHOD RMSE
% ================================================================

figure('Color','w');

bar(categorical(metricTable.Method),metricTable.RMSE);

ylabel('RMSE (BPM)');
xlabel('Estimator');
title('Heart-Rate Estimation RMSE by Method');

grid on;

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig2_Method_RMSE.png'));

%% ================================================================
% 26. FIGURE 3 - REFERENCE VS ADAPTIVE
% ================================================================

mask = validAdaptive;

figure('Color','w');

scatter( ...
    reference(mask), ...
    R.AdaptiveHR(mask), ...
    12, ...
    'filled');

hold on;

lims = [ ...
    min(reference(mask)), ...
    max(reference(mask))];

plot(lims,lims,'k--','LineWidth',1.5);

xlabel('Reference HR (BPM)');
ylabel('Adaptive estimated HR (BPM)');
title('Reference vs. Adaptive HR Estimation');

grid on;
axis equal;

xlim(lims);
ylim(lims);

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig3_Reference_vs_Adaptive.png'));

%% ================================================================
% 27. FIGURE 4 - ERROR DISTRIBUTION
% ================================================================

err = R.AdaptiveHR(mask)-reference(mask);

figure('Color','w');

histogram(err,30);

xlabel('Estimation Error (BPM)');
ylabel('Number of Records');
title('Adaptive HR Estimation Error Distribution');

grid on;

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig4_Error_Distribution.png'));

%% ================================================================
% 28. FIGURE 5 - RISK COVERAGE
% ================================================================

figure('Color','w');

methodsRC = unique(RC.Method,'stable');

for k = 1:numel(methodsRC)

    mm = RC.Method==methodsRC(k);

    plot( ...
        RC.Coverage(mm), ...
        RC.RiskMAE(mm), ...
        'LineWidth',1.5);

    hold on;

end

xlabel('Coverage');
ylabel('MAE (BPM)');
title('Risk-Coverage Curve');

legend(cellstr(methodsRC), ...
    'Location','best');

grid on;

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig5_Risk_Coverage.png'));

%% ================================================================
% 29. FIGURE 6 - RELIABILITY DISTRIBUTION
% ================================================================

relMask = success & isfinite(R.Reliability);

figure('Color','w');

histogram(R.Reliability(relMask),30);

xlabel('Reliability Score');
ylabel('Number of Records');
title('Distribution of Reliability Scores');

grid on;

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig6_Reliability.png'));

%% ================================================================
% 30. FIGURE 7 - MOTION VS ERROR
% ================================================================

mmask = validAdaptive & isfinite(R.MotionIndex);

figure('Color','w');

scatter( ...
    R.MotionIndex(mmask), ...
    abs(R.AdaptiveHR(mmask)-reference(mmask)), ...
    12, ...
    'filled');

xlabel('Motion Index');
ylabel('Absolute HR Error (BPM)');
title('Motion Contamination vs. Absolute Estimation Error');

grid on;

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig7_Motion_vs_Error.png'));

%% ================================================================
% 31. FIGURE 8 - BLAND ALTMAN
% ================================================================

figure('Color','w');

pred = R.AdaptiveHR(mask);
ref = reference(mask);

meanHR = (pred+ref)/2;
diffHR = pred-ref;

biasBA = mean(diffHR);
sdBA = std(diffHR);

loaUpper = biasBA + 1.96*sdBA;
loaLower = biasBA - 1.96*sdBA;

scatter(meanHR,diffHR,12,'filled');

hold on;

yline(biasBA,'k-','LineWidth',1.5);
yline(loaUpper,'k--','LineWidth',1.2);
yline(loaLower,'k--','LineWidth',1.2);

xlabel('Mean HR of Reference and Estimate (BPM)');
ylabel('Difference: Estimate - Reference (BPM)');
title('Bland-Altman Analysis');

grid on;

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig8_Bland_Altman.png'));

%% ================================================================
% 32. FIGURE 9 - QUALITY VS RELIABILITY
% ================================================================

qmask = success & ...
    isfinite(R.Quality) & ...
    isfinite(R.Reliability);

figure('Color','w');

boxchart( ...
    categorical(R.Quality(qmask)), ...
    R.Reliability(qmask));

xlabel('Dataset Quality Label');
ylabel('Reliability Score');
title('Reliability Score Across Reference Quality Classes');

grid on;

saveas(gcf, ...
    fullfile(FIG_ROOT,'Fig9_Quality_Reliability.png'));

%% ================================================================
% 33. FIGURE 10 - SUBJECT MAE
% ================================================================

if ~isempty(SubjectResults)

    figure('Color','w');

    histogram(SubjectResults.MAE,20);

    xlabel('Subject-Level MAE (BPM)');
    ylabel('Number of Subjects');
    title('Distribution of Subject-Level Adaptive MAE');

    grid on;

    saveas(gcf, ...
        fullfile(FIG_ROOT,'Fig10_Subject_MAE.png'));

end

%% ================================================================
% 34. FINAL SUMMARY
% ================================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' FINAL RESULTS SUMMARY\n');
fprintf('===============================================================\n');

fprintf('\n');

disp(metricTable);

fprintf('\nAdaptive method:\n');

adaptiveRow = metricTable( ...
    metricTable.Method=="Adaptive",:);

disp(adaptiveRow);

fprintf('\nResults directory:\n%s\n',RESULT_ROOT);

fprintf('\nFigures directory:\n%s\n',FIG_ROOT);

%% ================================================================
% 35. SAVE COMPLETE WORKSPACE
% ================================================================

save( ...
    fullfile(RESULT_ROOT,'BUT_PPG_JOURNAL_V2_COMPLETE.mat'), ...
    'R', ...
    'metricTable', ...
    'FoldResults', ...
    'SubjectResults', ...
    'BootstrapResults', ...
    'RC', ...
    'MC', ...
    'qualityResults', ...
    'hrResults', ...
    'motionResults', ...
    'ComparisonResults', ...
    'DATA_ROOT', ...
    'RESULT_ROOT', ...
    '-v7.3');

fprintf('\n');
fprintf('===============================================================\n');
fprintf(' PIPELINE COMPLETED SUCCESSFULLY\n');
fprintf('===============================================================\n');

fprintf('\nAll major outputs have been saved.\n');

%% ================================================================
% LOCAL FUNCTIONS
% ================================================================


function col = findColumn(T,candidates)

col = '';

vars = T.Properties.VariableNames;

for i = 1:numel(candidates)

    idx = strcmpi(vars,candidates{i});

    if any(idx)
        col = vars{find(idx,1)};
        return;
    end

end

end


%% ================================================================
% ROBUST LOCAL WFDB LOADER
% ================================================================

function [x,fs,ok,msg] = loadLocalWFDBRecord(recordBase)

x = [];
fs = [];
ok = false;
msg = '';

[parentDir,recordName,~] = fileparts(recordBase);

heaFile = [recordBase '.hea'];
datFile = [recordBase '.dat'];

%% Direct existence test

if ~exist(heaFile,'file')

    msg = sprintf( ...
        'Header file not found: %s',heaFile);

    return;

end

if ~exist(datFile,'file')

    msg = sprintf( ...
        'Data file not found: %s',datFile);

    return;

end

%% Store current working directory

oldDir = pwd;

cleanupObj = onCleanup(@() cd(oldDir)); %#ok<NASGU>

%% IMPORTANT:
% WFDB rdsamp can fail when supplied with a full Windows path
% depending on WFDB MATLAB configuration.
%
% Therefore temporarily move into the record directory and
% provide only the local record stem.

try

    cd(parentDir);

    [x,fs] = rdsamp(recordName);

    if isempty(x)

        msg = sprintf( ...
            'rdsamp returned empty signal for %s',recordName);

        return;

    end

    ok = true;

catch ME1

    %% Fallback: explicit current-directory path

    try

        [x,fs] = rdsamp(fullfile('.',recordName));

        if isempty(x)

            msg = sprintf( ...
                'Empty signal from fallback rdsamp for %s', ...
                recordName);

            return;

        end

        ok = true;

    catch ME2

        x = [];
        fs = [];

        msg = sprintf( ...
            ['WFDB loading failed for %s.\n' ...
             'Primary error: %s\n' ...
             'Fallback error: %s'], ...
             recordName,ME1.message,ME2.message);

    end

end

end


%% ================================================================
% PROCESS ONE RECORD
% ================================================================

function [out,ok,msg] = processOneRecord( ...
    DATA_ROOT, ...
    rid, ...
    FS_TARGET, ...
    LOW_HZ, ...
    HIGH_HZ, ...
    MIN_HR, ...
    MAX_HR, ...
    MIN_SECONDS, ...
    NFFT)

out = struct();

ok = false;

msg = '';

%% Record paths

recordDir = fullfile(DATA_ROOT,rid);

ppgBase = fullfile( ...
    recordDir,[rid '_PPG']);

accBase = fullfile( ...
    recordDir,[rid '_ACC']);

%% ---------------------------------------------------------------
% Load PPG
% ---------------------------------------------------------------

[x,fs,okPPG,msgPPG] = ...
    loadLocalWFDBRecord(ppgBase);

if ~okPPG

    msg = msgPPG;

    return;

end

%% ---------------------------------------------------------------
% Convert to single channel
% ---------------------------------------------------------------

if size(x,2)>1

    x = x(:,1);

end

x = double(x(:));

if isempty(fs) || ~isfinite(fs) || fs<=0

    msg = 'Invalid PPG sampling frequency.';

    return;

end

%% ---------------------------------------------------------------
% Remove NaN / Inf
% ---------------------------------------------------------------

x(~isfinite(x)) = 0;

%% ---------------------------------------------------------------
% Duration check
% ---------------------------------------------------------------

if numel(x)/fs < MIN_SECONDS

    msg = sprintf( ...
        'PPG shorter than MIN_SECONDS (%.2f s).', ...
        numel(x)/fs);

    return;

end

%% ---------------------------------------------------------------
% Resample
% ---------------------------------------------------------------

if abs(fs-FS_TARGET)>1e-6

    try

        x = resample(x,FS_TARGET,round(fs));

    catch

        tOld = (0:numel(x)-1)'/fs;

        tNew = (0:1/FS_TARGET:tOld(end))';

        x = interp1(tOld,x,tNew,'linear','extrap');

    end

end

fs = FS_TARGET;

%% ---------------------------------------------------------------
% Limit / detrend
% ---------------------------------------------------------------

x = x(:);

x = detrend(x,'constant');

%% Robust scaling

medX = median(x);

madX = median(abs(x-medX));

if madX>eps

    x = (x-medX)/(1.4826*madX);

else

    sdX = std(x);

    if sdX>eps
        x = (x-mean(x))/sdX;
    end

end

%% ---------------------------------------------------------------
% Bandpass filtering
% ---------------------------------------------------------------

try

    bp = designfilt( ...
        'bandpassiir', ...
        'FilterOrder',4, ...
        'HalfPowerFrequency1',LOW_HZ, ...
        'HalfPowerFrequency2',HIGH_HZ, ...
        'SampleRate',fs);

    xFilt = filtfilt(bp,x);

catch

    [b,a] = butter(3,[LOW_HZ HIGH_HZ]/(fs/2),'bandpass');

    xFilt = filtfilt(b,a,x);

end

xFilt = xFilt(:);

%% ---------------------------------------------------------------
% Load ACC if available
% ---------------------------------------------------------------

motionIndex = 0;

if exist([accBase '.hea'],'file') && ...
        exist([accBase '.dat'],'file')

    [acc,accFs,okACC,~] = ...
        loadLocalWFDBRecord(accBase);

    if okACC

        acc = double(acc);

        acc(~isfinite(acc)) = 0;

        motionIndex = calculateMotionIndex(acc,accFs);

    end

end

%% ---------------------------------------------------------------
% Candidate estimators
% ---------------------------------------------------------------

peakHR = estimatePeakHR( ...
    xFilt,fs,LOW_HZ,HIGH_HZ,MIN_HR,MAX_HR);

fftHR = estimateFFTHR( ...
    xFilt,fs,LOW_HZ,HIGH_HZ,MIN_HR,MAX_HR,NFFT);

psdHR = estimatePSDHR( ...
    xFilt,fs,LOW_HZ,HIGH_HZ,MIN_HR,MAX_HR);

acfHR = estimateACFHR( ...
    xFilt,fs,MIN_HR,MAX_HR);

harmonicHR = estimateHarmonicHR( ...
    xFilt,fs,LOW_HZ,HIGH_HZ,MIN_HR,MAX_HR,NFFT);

%% ---------------------------------------------------------------
% Quality features
% ---------------------------------------------------------------

periodicity = calculatePeriodicity(xFilt,fs);

spectralConcentration = ...
    calculateSpectralConcentration( ...
    xFilt,fs,LOW_HZ,HIGH_HZ,NFFT);

spectralEntropy = ...
    calculateSpectralEntropy( ...
    xFilt,fs,LOW_HZ,HIGH_HZ,NFFT);

candidateVector = [ ...
    peakHR ...
    fftHR ...
    psdHR ...
    acfHR ...
    harmonicHR];

finiteCandidates = candidateVector( ...
    isfinite(candidateVector));

if numel(finiteCandidates)>=2

    estimatorAgreement = ...
        1/(1+std(finiteCandidates));

else

    estimatorAgreement = 0;

end

%% ---------------------------------------------------------------
% Adaptive reliability
% ---------------------------------------------------------------

[adaptiveHR,reliability,selectedMethod] = ...
    adaptiveSelection( ...
    candidateVector, ...
    periodicity, ...
    spectralConcentration, ...
    spectralEntropy, ...
    estimatorAgreement, ...
    motionIndex);

%% ---------------------------------------------------------------
% Output
% ---------------------------------------------------------------

out.PeakHR = peakHR;
out.FFTHR = fftHR;
out.PSDHR = psdHR;
out.ACFHR = acfHR;
out.HarmonicHR = harmonicHR;
out.AdaptiveHR = adaptiveHR;

out.Periodicity = periodicity;
out.SpectralConcentration = spectralConcentration;
out.SpectralEntropy = spectralEntropy;
out.EstimatorAgreement = estimatorAgreement;
out.MotionIndex = motionIndex;

out.Reliability = reliability;

out.SelectedMethod = selectedMethod;

ok = true;

end


%% ================================================================
% PEAK ESTIMATOR
% ================================================================

function hr = estimatePeakHR(x,fs,lowHz,highHz,minHR,maxHR)

hr = NaN;

if isempty(x)
    return;
end

minDist = floor(fs*60/maxHR*0.55);

minProm = 0.10*std(x);

[pks,locs] = findpeaks( ...
    x, ...
    'MinPeakDistance',max(1,minDist), ...
    'MinPeakProminence',minProm);

if numel(locs)<2
    return;
end

ibi = diff(locs)/fs;

ibi = ibi(ibi>60/maxHR & ibi<60/minHR);

if isempty(ibi)
    return;
end

hr = 60/median(ibi);

if hr<minHR || hr>maxHR
    hr = NaN;
end

end


%% ================================================================
% FFT ESTIMATOR
% ================================================================

function hr = estimateFFTHR(x,fs,lowHz,highHz,minHR,maxHR,NFFT)

hr = NaN;

N = numel(x);

if N<10
    return;
end

x = x-mean(x);

w = hann(N);

X = abs(fft(x.*w,NFFT));

freq = (0:NFFT-1)'*fs/NFFT;

keep = freq>=lowHz & ...
       freq<=highHz & ...
       freq>=minHR/60 & ...
       freq<=maxHR/60;

if ~any(keep)
    return;
end

spec = X(keep);

freqK = freq(keep);

[~,idx] = max(spec);

hr = freqK(idx)*60;

end


%% ================================================================
% WELCH PSD
% ================================================================

function hr = estimatePSDHR(x,fs,lowHz,highHz,minHR,maxHR)

hr = NaN;

N = numel(x);

if N<32
    return;
end

segLen = min(N,round(8*fs));

if segLen<32
    return;
end

noverlap = floor(segLen/2);

nfft = max(1024,2^nextpow2(segLen));

try

    [pxx,f] = pwelch( ...
        x, ...
        hamming(segLen), ...
        noverlap, ...
        nfft, ...
        fs);

catch

    [pxx,f] = pwelch(x,[],[],nfft,fs);

end

keep = ...
    f>=lowHz & ...
    f<=highHz & ...
    f>=minHR/60 & ...
    f<=maxHR/60;

if ~any(keep)
    return;
end

[~,idx] = max(pxx(keep));

fk = f(keep);

hr = fk(idx)*60;

end


%% ================================================================
% ACF ESTIMATOR
% ================================================================

function hr = estimateACFHR(x,fs,minHR,maxHR)

hr = NaN;

x = x(:);

if numel(x)<20
    return;
end

x = x-mean(x);

maxLag = floor(fs*60/minHR);
minLag = ceil(fs*60/maxHR);

if maxLag>=numel(x)-1
    maxLag = numel(x)-1;
end

if minLag>=maxLag
    return;
end

acf = xcorr(x,maxLag,'coeff');

acf = acf(maxLag+1:end);

search = acf(minLag+1:maxLag+1);

if isempty(search)
    return;
end

[~,idx] = max(search);

lag = minLag+idx-1;

if lag<=0
    return;
end

hr = 60*fs/lag;

if hr<minHR || hr>maxHR
    hr = NaN;
end

end


%% ================================================================
% HARMONIC ESTIMATOR
% ================================================================

function hr = estimateHarmonicHR( ...
    x,fs,lowHz,highHz,minHR,maxHR,NFFT)

hr = NaN;

N = numel(x);

if N<20
    return;
end

x = x-mean(x);

w = hann(N);

X = abs(fft(x.*w,NFFT));

freq = (0:NFFT-1)'*fs/NFFT;

keep = ...
    freq>=lowHz & ...
    freq<=highHz;

freq = freq(keep);
spec = X(keep);

if isempty(spec)
    return;
end

%% Find local maxima

[pks,locs] = findpeaks(spec);

if isempty(pks)

    [~,idx] = max(spec);

    hr = freq(idx)*60;

    return;

end

candidateFreq = freq(locs);

candidatePower = pks;

%% Score candidate fundamental frequencies

scores = zeros(size(candidateFreq));

for i = 1:numel(candidateFreq)

    f0 = candidateFreq(i);

    if f0<minHR/60 || f0>maxHR/60
        continue;
    end

    score = 0;

    for h = 1:4

        target = h*f0;

        [d,idx] = min(abs(freq-target));

        if d <= 2*fs/NFFT

            score = score + ...
                spec(idx)/(h^0.7);

        end

    end

    scores(i) = score;

end

if all(scores<=0)

    [~,idx] = max(candidatePower);

    hr = candidateFreq(idx)*60;

else

    [~,idx] = max(scores);

    hr = candidateFreq(idx)*60;

end

if hr<minHR || hr>maxHR
    hr = NaN;
end

end


%% ================================================================
% PERIODICITY
% ================================================================

function p = calculatePeriodicity(x,fs)

p = 0;

x = x(:);

if numel(x)<20
    return;
end

x = x-mean(x);

maxLag = min( ...
    floor(2*fs), ...
    numel(x)-1);

minLag = max( ...
    ceil(fs*60/180), ...
    1);

if minLag>=maxLag
    return;
end

acf = xcorr(x,maxLag,'coeff');

acf = acf(maxLag+1:end);

search = acf(minLag+1:maxLag+1);

if isempty(search)
    return;
end

p = max(search);

p = max(0,min(1,p));

end


%% ================================================================
% SPECTRAL CONCENTRATION
% ================================================================

function c = calculateSpectralConcentration( ...
    x,fs,lowHz,highHz,NFFT)

c = 0;

N = numel(x);

if N<10
    return;
end

x = x-mean(x);

X = abs(fft(x.*hann(N),NFFT)).^2;

f = (0:NFFT-1)'*fs/NFFT;

keep = ...
    f>=lowHz & ...
    f<=highHz;

if ~any(keep)
    return;
end

s = X(keep);

if sum(s)<=eps
    return;
end

s = s/sum(s);

%% Use percentile energy concentration rather than a
% zero-padding-sensitive single-bin ratio.

sSorted = sort(s,'descend');

cumEnergy = cumsum(sSorted);

idx = find(cumEnergy>=0.50,1,'first');

if isempty(idx)
    idx = numel(sSorted);
end

c = 1 - (idx-1)/max(numel(sSorted)-1,1);

c = max(0,min(1,c));

end


%% ================================================================
% SPECTRAL ENTROPY
% ================================================================

function e = calculateSpectralEntropy( ...
    x,fs,lowHz,highHz,NFFT)

e = 1;

N = numel(x);

if N<10
    return;
end

x = x-mean(x);

P = abs(fft(x.*hann(N),NFFT)).^2;

f = (0:NFFT-1)'*fs/NFFT;

keep = ...
    f>=lowHz & ...
    f<=highHz;

if ~any(keep)
    return;
end

P = P(keep);

if sum(P)<=eps
    return;
end

P = P/sum(P);

P = P(P>0);

e = -sum(P.*log(P));

e = e/log(numel(P));

e = max(0,min(1,e));

end


%% ================================================================
% MOTION INDEX
% ================================================================

function motion = calculateMotionIndex(acc,fs)

motion = 0;

if isempty(acc)
    return;
end

acc = double(acc);

if size(acc,2)>3

    acc = acc(:,1:3);

end

acc(~isfinite(acc)) = 0;

%% Remove DC/gravity component using a moving average

win = max(3,round(fs*1.0));

for k = 1:size(acc,2)

    acc(:,k) = ...
        acc(:,k)-movmean(acc(:,k),win);

end

rmsMotion = sqrt(mean(acc.^2,1));

rmsMotion = rmsMotion(isfinite(rmsMotion));

if isempty(rmsMotion)
    return;
end

motionRaw = mean(rmsMotion);

%% Log compression prevents very large ACC amplitudes from
% dominating the reliability model.

motion = log1p(motionRaw);

%% Normalize to a practical [0,1] range

motion = motion/(1+motion);

motion = max(0,min(1,motion));

end


%% ================================================================
% ADAPTIVE SELECTION
% ================================================================

function [hr,reliability,selectedMethod] = ...
    adaptiveSelection( ...
    candidates, ...
    periodicity, ...
    spectralConcentration, ...
    spectralEntropy, ...
    estimatorAgreement, ...
    motion)

names = [ ...
    "Peak", ...
    "FFT", ...
    "PSD", ...
    "ACF", ...
    "Harmonic"];

valid = isfinite(candidates);

if ~any(valid)

    hr = NaN;
    reliability = 0;
    selectedMethod = "";

    return;

end

vals = candidates(valid);

%% Robust central estimate

central = median(vals);

%% Candidate agreement

if numel(vals)>=2

    deviations = abs(vals-central);

    candidateAgreement = ...
        exp(-median(deviations)/5);

else

    candidateAgreement = 0.5;

end

%% Individual candidate confidence

confidence = zeros(size(candidates));

for i = 1:numel(candidates)

    if ~isfinite(candidates(i))
        continue;
    end

    distance = abs(candidates(i)-central);

    confidence(i) = exp(-distance/6);

end

%% Physiological quality components

P = max(0,min(1,periodicity));

C = max(0,min(1,spectralConcentration));

E = max(0,min(1,1-spectralEntropy));

A = max(0,min(1, ...
    0.5*estimatorAgreement + ...
    0.5*candidateAgreement));

M = max(0,min(1,motion));

%% Reliability score

reliability = ...
    0.30*P + ...
    0.25*C + ...
    0.15*E + ...
    0.20*A + ...
    0.10*(1-M);

reliability = max(0,min(1,reliability));

%% Candidate score

% Fixed estimator priors are intentionally modest.
% They do not use the reference HR.

priors = [ ...
    1.00 ...
    0.95 ...
    1.00 ...
    0.90 ...
    0.95];

score = confidence.*priors;

%% Penalize extreme disagreement

score(~valid) = -inf;

[~,idx] = max(score);

if ~isfinite(score(idx))

    hr = central;

    selectedMethod = "Median";

    return;

end

selectedMethod = names(idx);

selected = candidates(idx);

%% Reliability-aware fusion:
% When the candidates are strongly clustered, use a weighted
% consensus instead of forcing one estimator.

if numel(vals)>=3 && ...
        candidateAgreement>0.75 && ...
        reliability>0.60

    w = confidence;

    w(~valid) = 0;

    if sum(w)>eps

        hr = sum(w.*candidates)/sum(w);

    else

        hr = selected;

    end

else

    hr = selected;

end

end


%% ================================================================
% METRICS
% ================================================================

function met = calculateMetrics(ref,pred)

ref = double(ref(:));
pred = double(pred(:));

valid = ...
    isfinite(ref) & ...
    isfinite(pred);

ref = ref(valid);
pred = pred(valid);

met = struct();

if isempty(ref)

    met.MAE = NaN;
    met.RMSE = NaN;
    met.MAPE = NaN;
    met.Bias = NaN;
    met.AER5 = NaN;
    met.AER10 = NaN;
    met.PearsonR = NaN;
    met.R2 = NaN;

    return;

end

err = pred-ref;

met.MAE = mean(abs(err));

met.RMSE = sqrt(mean(err.^2));

met.MAPE = ...
    mean(abs(err)./max(abs(ref),eps))*100;

met.Bias = mean(err);

met.AER5 = mean(abs(err)<=5)*100;

met.AER10 = mean(abs(err)<=10)*100;

if numel(ref)>=2 && ...
        std(ref)>eps && ...
        std(pred)>eps

    cc = corrcoef(ref,pred);

    met.PearsonR = cc(1,2);

else

    met.PearsonR = NaN;

end

SSres = sum((ref-pred).^2);

SStot = sum((ref-mean(ref)).^2);

if SStot>eps

    met.R2 = 1-SSres/SStot;

else

    met.R2 = NaN;

end

end


%% ================================================================
% SUBJECT-LEVEL BOOTSTRAP MAE
% ================================================================

function [bootMean,CI] = ...
    subjectBootstrapMAE( ...
    subjectID,ref,pred,NBOOT,seed)

rng(seed);

subjectID = subjectID(:);
ref = ref(:);
pred = pred(:);

subjects = unique(subjectID);

S = numel(subjects);

if S==0

    bootMean = NaN;
    CI = [NaN NaN];

    return;

end

subjectMAE = nan(S,1);

for s = 1:S

    mask = subjectID==subjects(s);

    subjectMAE(s) = ...
        mean(abs(pred(mask)-ref(mask)));

end

bootMean = mean(subjectMAE);

B = nan(NBOOT,1);

for b = 1:NBOOT

    idx = randi(S,S,1);

    B(b) = mean(subjectMAE(idx));

end

CI = prctile(B,[2.5 97.5]);

end


%% ================================================================
% RISK-COVERAGE TABLE
% ================================================================

function RC = riskCoverageTable(R)

methods = { ...
    'Peak', ...
    'FFT', ...
    'PSD', ...
    'ACF', ...
    'Harmonic', ...
    'Adaptive'};

RC = table();

for m = 1:numel(methods)

    method = methods{m};

    switch method

        case 'Peak'
            pred = R.PeakHR;

        case 'FFT'
            pred = R.FFTHR;

        case 'PSD'
            pred = R.PSDHR;

        case 'ACF'
            pred = R.ACFHR;

        case 'Harmonic'
            pred = R.HarmonicHR;

        case 'Adaptive'
            pred = R.AdaptiveHR;

    end

    valid = ...
        R.Status=="SUCCESS" & ...
        isfinite(pred) & ...
        isfinite(R.ReferenceHR) & ...
        isfinite(R.Reliability);

    if sum(valid)<2
        continue;
    end

    p = pred(valid);
    r = R.ReferenceHR(valid);
    rel = R.Reliability(valid);

    %% Highest reliability first

    [~,order] = sort(rel,'descend');

    p = p(order);
    r = r(order);

    n = numel(p);

    coverage = (1:n)'/n;

    risk = nan(n,1);

    absErr = abs(p-r);

    for k = 1:n

        risk(k) = mean(absErr(1:k));

    end

    methodCol = repmat(string(method),n,1);

    temp = table( ...
        methodCol, ...
        coverage, ...
        risk, ...
        'VariableNames', ...
        {'Method','Coverage','RiskMAE'});

    RC = [RC;temp];

end

end


%% ================================================================
% MATCHED COVERAGE
% ================================================================

function MC = matchedCoverageTable(R,targetCoverage)

methods = { ...
    'Peak', ...
    'FFT', ...
    'PSD', ...
    'ACF', ...
    'Harmonic', ...
    'Adaptive'};

allValid = ...
    R.Status=="SUCCESS" & ...
    isfinite(R.ReferenceHR) & ...
    isfinite(R.Reliability);

predMatrix = [ ...
    R.PeakHR ...
    R.FFTHR ...
    R.PSDHR ...
    R.ACFHR ...
    R.HarmonicHR ...
    R.AdaptiveHR];

allValid = ...
    allValid & ...
    all(isfinite(predMatrix),2);

idx = find(allValid);

if isempty(idx)

    MC = table();

    return;

end

%% Rank the common records using proposed reliability

[~,order] = sort( ...
    R.Reliability(idx), ...
    'descend');

idx = idx(order);

nSelect = max(1,round(targetCoverage*numel(idx)));

idx = idx(1:nSelect);

MC = table();

for m = 1:numel(methods)

    pred = predMatrix(idx,m);

    ref = R.ReferenceHR(idx);

    met = calculateMetrics(ref,pred);

    row = table( ...
        string(methods{m}), ...
        numel(idx), ...
        numel(idx)/numel(allValid), ...
        met.MAE, ...
        met.RMSE, ...
        met.MAPE, ...
        met.Bias, ...
        met.AER5, ...
        met.AER10, ...
        met.PearsonR, ...
        met.R2, ...
        'VariableNames',{ ...
        'Method','N','Coverage','MAE','RMSE','MAPE','Bias', ...
        'AER5','AER10','PearsonR','R2'});

    MC = [MC;row];

end

end


%% ================================================================
% STRATIFIED METRICS
% ================================================================

function T = stratifiedMetrics( ...
    variable, ...
    bins, ...
    labels, ...
    pred, ...
    ref, ...
    valid)

T = table();

variable = double(variable(:));
pred = double(pred(:));
ref = double(ref(:));
valid = logical(valid(:));

for b = 1:numel(labels)

    mask = valid & ...
        variable>=bins(b) & ...
        variable<bins(b+1);

    if sum(mask)<1
        continue;
    end

    met = calculateMetrics(ref(mask),pred(mask));

    row = table( ...
        string(labels{b}), ...
        sum(mask), ...
        met.MAE, ...
        met.RMSE, ...
        met.MAPE, ...
        met.Bias, ...
        met.AER5, ...
        met.AER10, ...
        met.PearsonR, ...
        met.R2, ...
        'VariableNames',{ ...
        'Stratum','N','MAE','RMSE','MAPE','Bias', ...
        'AER5','AER10','PearsonR','R2'});

    T = [T;row];

end

end


%% ================================================================
% PAIRED COMPARISON TABLE
% ================================================================

function T = pairedComparisonTable(R)

methods = { ...
    'Peak', ...
    'FFT', ...
    'PSD', ...
    'ACF', ...
    'Harmonic', ...
    'Adaptive'};

preds = { ...
    R.PeakHR, ...
    R.FFTHR, ...
    R.PSDHR, ...
    R.ACFHR, ...
    R.HarmonicHR, ...
    R.AdaptiveHR};

T = table();

%% Compare each fixed estimator against Adaptive

for m = 1:5

    p1 = preds{m};
    p2 = preds{6};

    mask = ...
        R.Status=="SUCCESS" & ...
        isfinite(R.ReferenceHR) & ...
        isfinite(p1) & ...
        isfinite(p2);

    if sum(mask)<5
        continue;
    end

    ref = R.ReferenceHR(mask);

    e1 = abs(p1(mask)-ref);
    e2 = abs(p2(mask)-ref);

    d = e1-e2;

    medianDifference = median(d);

    %% Sign-rank style approximation using signed ranks

    pValue = approximateWilcoxon(d);

    effect = ...
        medianDifference/max(median(abs(d)),eps);

    row = table( ...
        string(methods{m}), ...
        "Adaptive", ...
        sum(mask), ...
        medianDifference, ...
        pValue, ...
        effect, ...
        'VariableNames',{ ...
        'Method1','Method2','N', ...
        'MedianAbsoluteErrorDifference', ...
        'PValue','EffectSize'});

    T = [T;row];

end

end


%% ================================================================
% APPROXIMATE WILCOXON SIGNED-RANK TEST
% ================================================================

function p = approximateWilcoxon(d)

d = d(:);

d = d(isfinite(d) & d~=0);

n = numel(d);

if n<5

    p = NaN;

    return;

end

absD = abs(d);

[sortedVals,order] = sort(absD);

ranks = zeros(n,1);

i = 1;

while i<=n

    j = i;

    while j<n && sortedVals(j+1)==sortedVals(i)
        j = j+1;
    end

    ranks(order(i:j)) = mean(i:j);

    i = j+1;

end

Wplus = sum(ranks(d>0));

Wminus = sum(ranks(d<0));

W = min(Wplus,Wminus);

mu = n*(n+1)/4;

sigma = sqrt(n*(n+1)*(2*n+1)/24);

if sigma<=0

    p = NaN;

    return;

end

z = (W-mu)/sigma;

p = erfc(abs(z)/sqrt(2));

end