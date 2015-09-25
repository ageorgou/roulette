% This file shows how to invoke one of the samplers, starting with the
% specification of the system, the observations and the sampler parameters.
% To use different samplers (other than the Gibbs sampler with truncation,
% used here) or with other models, adjust this file or consult the other
% examples --- the setup is similar.

%% Define the system (via its reactions) and the parameters:
addpath('../util'); % for helper functions
addpath('../sampler'); % for the samplers

obs = load('../obsPredPreyRT_true'); % observation file

% First, the jump or update vector for each reaction:
updates = [+1 0; -1 0; 0 +1; 0 -1];

% Then the reaction propensities (the rate of reaction i is k_i * f_i(x),
% where x is the state, k_i is the parameter for reaction i and f_i is its
% propensity..
% Ensure that the functions f_i are vectorised, i.e. that they can be
% called with multiple states at the same time and return the rate from
% each of those states. This is to accelerate the computation of generator
% matrices.
rates = { @(x) x(:,1) .* x(:,2); ... % birth of predator
    @(x) x(:,1); ... % death of predator
    @(x) x(:,2); ... % birth of prey
    @(x) x(:,1) .* x(:,2)}; % death of prey

% Collect the updates and rates in a structure, and add the priors for each
% parameter, i.e. the parameters of its Gamma distribution.
% Note that the code uses the parameterisation of a Gamma in terms of shape
% and rate parameters, while MATLAB uses shape and scale (1/rate).
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


%% Run the sampler

% Number of samples to take. The samplers print a diagnostic message every
% PRINT_EVERY samples.
nSamples = 5000;
global PRINT_EVERY
PRINT_EVERY = 500;

fprintf('Starting modified Gibbs sampler...\n');
tic;
samplesG = sampler(nSamples,obs,reactions); % Gibbs sampler with truncation
timeElapsed = toc;
fprintf('Time elapsed: %f s\n',timeElapsed);

%% Write samples and time to file

fid = fopen('experiments/predPreyGibbs/samplesPredPreyGibbs','w');
formatString = ['%f' repmat(' %f',1,length(reactions)-1) '\n'];
fprintf(fid,formatString,samplesG');
fclose(fid);

fid = fopen('experiments/predPreyGibbs/timePredPreyGibbs','w');
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