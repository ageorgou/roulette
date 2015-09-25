%% MH sampler for params (predator-prety model)
addpath('../util'); % for helper functions
addpath('../expmv/'); % for efficient matrix exponentiation
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

nSamples = 5000; % Total samples to take
% A few word on the globals:
% The proposal distribution is Gaussian with standard deviation PROP_STD(i)
% in dimension i.
% A diagnostic is printed every PRINT_EVERY samples.
% ROUL_CALLS and ROUL_TERMS keep track of the total number of calls to the
% Russian Roulette algorithm and the total terms taken.
global PROP_STD PRINT_EVERY ROUL_CALLS ROUL_TERMS;
PROP_STD = [0.00005 0.00005 0.00005 0.00005];
PRINT_EVERY = 200;
ROUL_CALLS = 0;
ROUL_TERMS = 0;

fprintf('Starting M-H sampler...\n');
tic;
samplesMH = MHSampler(nSamples,obs,reactions);
timeElapsed = toc;
fprintf('Time elapsed: %f s\n',timeElapsed);
fprintf('Average terms per roulette: %f\n', ROUL_TERMS / ROUL_CALLS);

%% Write samples and time to file

fid = fopen('experiments/predPreyMH_RT/samplesPredPreyMH','w');
formatString = ['%f' repmat(' %f',1,length(reactions)-1) '\n'];
fprintf(fid,formatString,samplesMH');
fclose(fid);

fid = fopen('experiments/predPreyMH_RT/timePredPreyMH','w');
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
% h = plotSamplesWithPriors(samplesMH,reactions,ranges,trueVals,parNames);
% for ii = 1:length(reactions)
%     print(h(ii),['histPar' num2str(ii) '.eps'],'-depsc');
% end