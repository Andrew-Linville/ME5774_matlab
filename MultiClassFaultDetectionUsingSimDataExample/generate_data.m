% 1. Extract supporting files and load parameters
unzip('pdmRecipPump_supportingfiles.zip');
pdmRecipPump_Parameters;
CAT_Pump_1051_DataFile_imported;

% 2. Define fault parameter variations
numParValues = 10;
leak_area_set_factor = linspace(0.00,0.036,numParValues);
leak_area_set = leak_area_set_factor*TRP_Par.Check_Valve.In.Max_Area;
leak_area_set = max(leak_area_set,1e-9); 
blockinfactor_set = linspace(0.8,0.53,numParValues);
bearingfactor_set = linspace(0,6e-4,numParValues);

% 3. Create simulation input combinations
nPerGroup = 100;
rng('default');

% Generate fault combinations (No fault, Single, Double, Triple)
leakArea = repmat(leak_area_set(1),nPerGroup,1);
blockingFactor = repmat(blockinfactor_set(1),nPerGroup,1);
bearingFactor = repmat(bearingfactor_set(1),nPerGroup,1);

% Single faults
idx = ceil(10*rand(nPerGroup,1));
leakArea = [leakArea; leak_area_set(idx)'];
blockingFactor = [blockingFactor;repmat(blockinfactor_set(1),nPerGroup,1)];
bearingFactor = [bearingFactor;repmat(bearingfactor_set(1),nPerGroup,1)];
idx = ceil(10*rand(nPerGroup,1));
leakArea = [leakArea; repmat(leak_area_set(1),nPerGroup,1)];
blockingFactor = [blockingFactor;blockinfactor_set(idx)'];
bearingFactor = [bearingFactor;repmat(bearingfactor_set(1),nPerGroup,1)];
idx = ceil(10*rand(nPerGroup,1));
leakArea = [leakArea; repmat(leak_area_set(1),nPerGroup,1)];
blockingFactor = [blockingFactor;repmat(blockinfactor_set(1),nPerGroup,1)];
bearingFactor = [bearingFactor;bearingfactor_set(idx)'];

% Double faults
idxA = ceil(10*rand(nPerGroup,1)); idxB = ceil(10*rand(nPerGroup,1));
leakArea = [leakArea; leak_area_set(idxA)'];
blockingFactor = [blockingFactor;blockinfactor_set(idxB)'];
bearingFactor = [bearingFactor;repmat(bearingfactor_set(1),nPerGroup,1)];
idxA = ceil(10*rand(nPerGroup,1)); idxB = ceil(10*rand(nPerGroup,1));
leakArea = [leakArea; leak_area_set(idxA)'];
blockingFactor = [blockingFactor;repmat(blockinfactor_set(1),nPerGroup,1)];
bearingFactor = [bearingFactor;bearingfactor_set(idxB)'];
idxA = ceil(10*rand(nPerGroup,1)); idxB = ceil(10*rand(nPerGroup,1));
leakArea = [leakArea; repmat(leak_area_set(1),nPerGroup,1)];
blockingFactor = [blockingFactor;blockinfactor_set(idxA)'];
bearingFactor = [bearingFactor;bearingfactor_set(idxB)'];

% Triple faults
idxA = ceil(10*rand(nPerGroup,1)); idxB = ceil(10*rand(nPerGroup,1)); idxC = ceil(10*rand(nPerGroup,1));
leakArea = [leakArea; leak_area_set(idxA)'];
blockingFactor = [blockingFactor;blockinfactor_set(idxB)'];
bearingFactor = [bearingFactor;bearingfactor_set(idxC)'];

% 4. Build SimulationInput objects
mdl = 'pdmRecipPump';
load_system(mdl); % Load the model into memory silently to prevent the block diagram error

for ct = numel(leakArea):-1:1
    simInput(ct) = Simulink.SimulationInput(mdl);
    simInput(ct) = setVariable(simInput(ct),'leak_cyl_area_WKSP',leakArea(ct));
    simInput(ct) = setVariable(simInput(ct),'block_in_factor_WKSP',blockingFactor(ct));
    simInput(ct) = setVariable(simInput(ct),'bearing_fault_frict_WKSP',bearingFactor(ct));
    simInput(ct) = setVariable(simInput(ct),'noise_seed_offset_WKSP',ct-1);
end

% 5. Run Ensemble
if isfolder('./Data')
    delete('./Data/*.mat')
end
[ok,e] = generateSimulationEnsemble(simInput,fullfile('.','Data'),'UseParallel',false,'ShowProgress',false);
ens = simulationEnsembleDatastore(fullfile('.','Data'));

% 6. Extract Features
ens.DataVariables = ["qOut_meas"; "SimulationInput"];
ens.ConditionVariables = ["LeakFault","BlockingFault","BearingFault"];
reset(ens);

while hasdata(ens)
    data = read(ens);
    
    % Preprocess
    tMin = seconds(0.8);
    flow = data.qOut_meas{1};
    flow = flow(flow.Time >= tMin,:);
    flow.Time = flow.Time - flow.Time(1);
    flow = retime(flow,'regular','linear','TimeStep',seconds(1e-3));
    fA = flow;
    fA.Data = fA.Data - mean(fA.Data);
    [flowP,flowF] = pspectrum(fA,'FrequencyLimits',[2 250]);
    
    % Get Fault Values
    simin = data.SimulationInput{1};
    vars = {simin.Variables.Name};
    LeakFault = simin.Variables(strcmp(vars,'leak_cyl_area_WKSP')).Value;
    BlockingFault = simin.Variables(strcmp(vars,'block_in_factor_WKSP')).Value;
    BearingFault = simin.Variables(strcmp(vars,'bearing_fault_frict_WKSP')).Value;
    faultValues = {'LeakFault', LeakFault, 'BlockingFault', BlockingFault, 'BearingFault', BearingFault};
    
    % Extract CI
    pMax = max(flowP); fPeak = flowF(flowP==pMax);
    pLow = sum(flowP(flowF >= 10 & flowF <= 20));
    pMid = sum(flowP(flowF >= 40 & flowF <= 60));
    pHigh = sum(flowP(flowF >= 100));
    [pKur,fKur] = pkurtosis(flow); pKur = fKur(pKur == max(pKur));
    csFlowRange = max(cumsum(flow.Data))-min(cumsum(flow.Data));
    
    feat = {'qMean', mean(flow.Data), 'qVar', var(flow.Data), 'qSkewness', skewness(flow.Data), ...
            'qKurtosis', kurtosis(flow.Data), 'qPeak2Peak', peak2peak(flow.Data), ...
            'qCrest', peak2rms(flow.Data), 'qRMS', rms(flow.Data), 'qMAD', mad(flow.Data), ...
            'qCSRange',csFlowRange, 'fPeak', fPeak, 'pLow', pLow, 'pMid', pMid, 'pHigh', pHigh, 'pKurtosis', pKur(1)};
            
    writeToLastMemberRead(ens, [faultValues, feat]);
end

% 7. Format and Export to CSV
reset(ens);
ens.SelectedVariables = ["fPeak","pLow","pMid","pHigh","pKurtosis",...
    "qMean","qVar","qSkewness","qKurtosis","qPeak2Peak","qCrest",...
    "qRMS","qMAD","qCSRange","LeakFault","BlockingFault","BearingFault"];
data = gather(tall(ens));

data.LeakFlag = data.LeakFault > 1e-6;
data.BlockingFlag = data.BlockingFault < 0.8;
data.BearingFlag = data.BearingFault > 0;
data.CombinedFlag = data.LeakFlag + 2*data.BlockingFlag + 4*data.BearingFlag;

data(:, {'LeakFault', 'BlockingFault', 'BearingFault', 'LeakFlag', 'BlockingFlag', 'BearingFlag'}) = [];
writetable(data, 'pump_features.csv');