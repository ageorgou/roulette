%% Simple Gibbs sampler for params (predator-prety model, imposed truncation)
addpath('../util'); % for helper functions
addpath('../sampler'); % for the samplers

obs = load('../obsPredPreyRT_true');

updates = [+1 0; -1 0; 0 +1; 0 -1];
rates = { @(x) x(:,1) .* x(:,2); ... % birth of predator
    @(x) x(:,1); ... % death of predator
    @(x) x(:,2); ... % birth of prey
    @(x) x(:,1) .* x(:,2)}; % death of prey

reactions = struct('name',{'predB' 'predD','preyB','preyD'}, ...
    'propensity',[], 'update',[],'aPrior',[],'bPrior',[]);
reactions(1).propensity = rates{1};
reactions(1).update = updates(1,:);
reactions(1).aPrior = 4;
reactions(1).bPrior = 10000;
reactions(2).propensity = rates{2};
reactions(2).update = updates(2,:);
reactions(2).aPrior = 4;
reactions(2).bPrior = 10000;
reactions(3).propensity = rates{3};
reactions(3).update = updates(3,:);
reactions(3).aPrior = 4;
reactions(3).bPrior = 10000;
reactions(4).propensity = rates{4};
reactions(4).update = updates(4,:);
reactions(4).aPrior = 4;
reactions(4).bPrior = 10000;

% Build the state-space, allowing up to 100 of either species:
upperLimit = 100;
x = repmat((0:upperLimit),upperLimit+1,1);
y = reshape(repmat((0:upperLimit)',1,upperLimit+1),upperLimit+1,upperLimit+1);
nStates = (upperLimit+1) * (upperLimit+1);
space = [(1:nStates)' x(:) y(:)];

nSamples = 5000;
global PRINT_EVERY
PRINT_EVERY = 100;

fprintf('Starting Rao-Teh sampler (with imposed truncation)...\n');
tt = tic;
samplesRT = parameterRTSampler(nSamples,obs,reactions,space);
timeElapsed = toc(tt);
fprintf('Time elapsed: %f s\n',timeElapsed);

%% Write samples and time to file

fid = fopen('experiments/predPreyRT/samplesPredPreyRT','w');
formatString = ['%f' repmat(' %f',1,length(reactions)-1) '\n'];
fprintf(fid,formatString,samplesRT');
fclose(fid);

fid = fopen('experiments/predPreyRT/timePredPreyGibbsRT','w');
fprintf(fid,'Time elapsed: %f s\n',timeElapsed);
fclose(fid);

fprintf('Done.\n');

% Quit MATLAB, but only if called from the command line
if ~usejava('desktop')
    fprintf('Exiting...\n');
    exit;
end

%% Plot results
% trueVals = [0.0001 0.0005 0.0005 0.0001];
% ranges = [0 0.001; 0 0.001; 0 001; 0 0.001];
% parNames = {'Birth rate of predator','Death rate of predator', ...
%     'Birth rate of prey','Death rate of prey'};
% h = plotSamplesWithPriors(samplesG,reactions,ranges,trueVals,parNames);
% for ii = 1:length(reactions)
%     print(h(ii),['histPar' num2str(ii) '.eps'],'-depsc');
% end
