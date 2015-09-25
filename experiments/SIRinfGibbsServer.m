%% Gibbs sampler for params with truncation (infinite SIR model)
addpath('../util'); % for helper functions
addpath('../sampler'); % for the samplers

obs = load('../obsSIRinf');

updates = [-1 +1 0; 0 -1 +1; +1 0 0];

% parameter values used for data: [0.4 0.5 0.4];

reactions = struct('name',{'birth' 'death'},'propensity',[], ...
    'update',[],'aPrior',[],'bPrior',[]);
reactions(1).propensity = @(x) x(:,1) .* x(:,2);
reactions(1).update = updates(1,:);
reactions(1).aPrior = 1.5;
reactions(1).bPrior = 5;

reactions(2).propensity = @(x) x(:,2);
reactions(2).update = updates(2,:);
reactions(2).aPrior = 1.5;
reactions(2).bPrior = 5;

reactions(3).propensity = @(x) ones(size(x,1),1);
reactions(3).update = updates(3,:);
reactions(3).aPrior = 1.5;
reactions(3).bPrior = 5;

nSamples = 5000;
global PRINT_EVERY
PRINT_EVERY = 200;

fprintf('Starting modified Gibbs sampler...\n');
tic;
samplesG = sampler(nSamples,obs,reactions);
timeElapsed = toc;
fprintf('Time elapsed: %f s\n',timeElapsed);

%% Write samples and time to file
fid = fopen('experiments/SIRinfGibbs/samplesSIRinfGibbs','w');
formatString = ['%f' repmat(' %f',1,length(reactions)-1) '\n'];
fprintf(fid,formatString,samplesG');
fclose(fid);

fid = fopen('experiments/SIRinfGibbs/timeSIRinfGibbs','w');
fprintf(fid,'Time elapsed: %f s\n',timeElapsed);
fclose(fid);

fprintf('Done.\n');

% Quit MATLAB, but only if called from the command line
if ~usejava('desktop')
    fprintf('Exiting...\n');
    exit;
end
