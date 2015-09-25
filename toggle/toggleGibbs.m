%% Gibbs sampler for params
addpath('../util'); % for helper functions
obs = load('../obsTogglePartial');

%trueValues = [5 5 0.25 0.25 0.1 0.1 0.1 0.1];
nReacts = 8;

rExp = 0.1; % in the exponents of the deactivation rates, assumed known

% Species order: Repr1 Repr2 Prom1 Prom2

% Prior hyperparameters used for first experiment in comments
reactions = struct('name',[],'propensity',[], ...
    'update',[],'aPrior',[],'bPrior',[]);

% Protein 1 expression
reactions(1).name = 'repr1Expr';
reactions(1).propensity = @(x) x(:,3);
reactions(1).update = [1 0 0 0];
reactions(1).aPrior = 1.5; %peak near 1.5 
reactions(1).bPrior = 1/3;

% Protein 2 expression
reactions(2).name = 'repr2Expr';
reactions(2).propensity = @(x) x(:,4);
reactions(2).update = [0 1 0 0];
reactions(2).aPrior = 1.5;
reactions(2).bPrior = 1/3;

% Protein 1 degradation
reactions(3).name = 'repr1Degr';
reactions(3).propensity = @(x) x(:,1);
reactions(3).update = [-1 0 0 0];
reactions(3).aPrior = 3; %5;
reactions(3).bPrior = 1/0.2; %1/0.1;

% Protein 2 degradation
reactions(4).name = 'repr2Degr';
reactions(4).propensity = @(x) x(:,2);
reactions(4).update = [0 -1 0 0];
reactions(4).aPrior = 3; %5;
reactions(4).bPrior = 1/0.2; %1/0.1;

% Gene 1 activation
reactions(5).name = 'prom1Act';
reactions(5).propensity = @(x) ones(size(x,1),1) - x(:,3);
reactions(5).update = [0 0 1 0];
reactions(5).aPrior = 3; %2.3;
reactions(5).bPrior = 1/0.2; %1/0.2;

% Gene 2 activation
reactions(6).name = 'prom2Act';
reactions(6).propensity = @(x) ones(size(x,1),1) - x(:,4);
reactions(6).update = [0 0 0 1];
reactions(6).aPrior = 3; %2.3
reactions(6).bPrior = 1/0.2; %1/0.2;

% Gene 1 deactivation
reactions(7).name = 'prom1Deact';
reactions(7).propensity = @(x) x(:,3) .* exp(rExp*x(:,2));
reactions(7).update = [0 0 -1 0];
reactions(7).aPrior = 3; %1.2;
reactions(7).bPrior = 1/0.2; %1/0.1;

% Gene 2 deactivation
reactions(8).name = 'prom2Deact';
reactions(8).propensity = @(x) x(:,4) .* exp(rExp*x(:,1));
reactions(8).update = [0 0 0 -1];
reactions(8).aPrior = 3; %1.2;
reactions(8).bPrior = 1/0.2; %1/0.1;

nSamples = 5000;
global PRINT_EVERY SAVE_EVERY;
PRINT_EVERY = 200;
SAVE_EVERY = 500;

folder = 'experiments/toggleGibbs2'; % second experiment

fprintf('Starting modified Gibbs sampler...\n');
tic;
samplesG = sampler_toggle(nSamples,obs,reactions,folder);
timeElapsed = toc;
fprintf('Time elapsed: %f s\n',timeElapsed);

%% Write samples and time to file
fid = fopen([folder '/samplesToggleGibbs'],'w');
formatString = ['%f' repmat(' %f',1,length(reactions)-1) '\n'];
fprintf(fid,formatString,samplesG');
fclose(fid);

fid = fopen([folder '/timeToggleGibbs'],'w');
fprintf(fid,'Time elapsed: %f s\n',timeElapsed);
fclose(fid);

fprintf('Done.\n');

% Quit MATLAB, but only if called from the command line
if ~usejava('desktop')
    fprintf('Exiting...\n');
    exit;
end

%% Plot results
% trueVals = [0.02 0.05 1 0.05];
% ranges = [0 0.1; 0 0.1; 0 5; 0 0.1];
% parNames = {'Birth rate of predator','Death rate of predator', ...
%     'Birth rate of prey','Death rate of prey'};
% h = plotSamplesWithPriors(samplesG,reactions,ranges,trueVals,parNames);
% for ii = 1:length(reactions)
%     print(h(ii),['histPar' num2str(ii) '.eps'],'-depsc');
% end
