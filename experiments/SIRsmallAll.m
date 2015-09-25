%% Set up for the finite SIR model
addpath('../util'); % for helper functions
addpath('../expmv/'); % for efficient matrix exponentiation
addpath('../sampler'); % for the samplers

obs = load('../obsSIR');

updates = [-1 +1 0; 0 -1 +1; +1 0 0];
%trueValues = [0.4 0.5];

reactions = struct('name',{'spread' 'recov'},'propensity',[], ...
    'update',[],'aPrior',[],'bPrior',[]);

reactions(1).propensity = @(x) x(:,1) .* x(:,2);
reactions(1).update = updates(1,:);
reactions(1).aPrior = 1.5;
reactions(1).bPrior = 5;

reactions(2).propensity = @(x) x(:,2);
reactions(2).update = updates(2,:);
reactions(2).aPrior = 1.5;
reactions(2).bPrior = 5;

nSamples = 5000;
global PRINT_EVERY
PRINT_EVERY = 500;

%% Rao-Teh sampler
% space = makeSIRStatespace([10 5 0]);
% 
% t0 = tic;
% samplesRT = parameterRTSampler(nSamples,obs,reactions,space);
% timeElapsed = toc(t0);
% fprintf('Time elapsed: %f s\n',timeElapsed);
% 
% fid = fopen('experiments/SIRsmallRT/samplesSIRsmallRT','w');
% formatString = ['%f' repmat(' %f',1,length(reactions)-1) '\n'];
% fprintf(fid,formatString,samplesRT');
% fclose(fid);
% 
% fid = fopen('experiments/SIRsmallRT/timeSIRsmallRT','w');
% fprintf(fid,'Time elapsed: %f s\n',timeElapsed);
% fclose(fid);
% fprintf('Done with Rao-Teh parameter sampler.\n');
% 
%% Metropolis-Hastings with Roulette
global PROP_STD ROUL_CALLS ROUL_TERMS;
PROP_STD = [0.2 0.2];

t0 = tic;
samplesMH = MHSampler(nSamples,obs,reactions);
timeElapsed = toc(t0);
fprintf('Time elapsed: %f s\n',timeElapsed);
fprintf('Average terms per roulette: %f\n', ROUL_TERMS / ROUL_CALLS);

fid = fopen('experiments/SIRsmallMH/samplesSIRsmallMH','w');
formatString = ['%f' repmat(' %f',1,length(reactions)-1) '\n'];
fprintf(fid,formatString,samplesMH');
fclose(fid);

fid = fopen('experiments/SIRsmallMH/timeSIRsmallMH','w');
fprintf(fid,'Time elapsed: %f s\n',timeElapsed);
fclose(fid);
fprintf('Done with Metropolis-Hastings Roulette sampler.\n');

%% Gibbs with Roulette

% t0 = tic;
% samplesG = sampler(nSamples,obs,reactions);
% timeElapsed = toc(t0);
% fprintf('Time elapsed: %f s\n',timeElapsed);
% 
% fid = fopen('experiments/SIRsmallGibbs/samplesSIRsmallGibbs','w');
% formatString = ['%f' repmat(' %f',1,length(reactions)-1) '\n'];
% fprintf(fid,formatString,samplesG');
% fclose(fid);
% 
% fid = fopen('experiments/SIRsmallGibbs/timeSIRsmallGibbs','w');
% fprintf(fid,'Time elapsed: %f s\n',timeElapsed);
% fclose(fid);
% fprintf('Done with Gibbs Roulette sampler.\n');

% Quit MATLAB, but only if called from the command line
if ~usejava('desktop')
    fprintf('Exiting...\n');
    exit;
end
