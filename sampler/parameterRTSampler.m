function parSamples = parameterRTSampler(nSamples,obs,reactions,states,outputFolder)
% PARAMETERRTSAMPLER Gibbs sampler for models with finite state-spaces.
%
% Following "Fast MCMC Sampling for Markov Jump Processes and Extensions"
% (Vinayak Rao and Yee Whye Teh, JMLR 2013)

global PRINT_EVERY SAVE_EVERY;

if isempty(SAVE_EVERY) %if it has not been defined, avoid saving
    SAVE_EVERY = nSamples + 1; %quick workaround
else
    if nargin < 4 || isempty(outputFolder)
        outputFolder = './';
    elseif outputFolder(end) ~= '/'
        outputFolder = [outputFolder '/'];
    end
end

R = length(reactions); % number of reactions / parameters
parSamples = zeros(nSamples,R);

updates = cell2mat({reactions.update}');
rateFuncs = {reactions.propensity};
aPrior = [reactions.aPrior];
bPrior = [reactions.bPrior];

aPost = aPrior; bPost = bPrior;

n = 1;
nResampled = 0;
while n <= nSamples
    params = sampleGamma(aPost,bPost);
    % Make generator
    A = makeGeneratorRest(states,updates,rateFuncs,params);
    % Draw a new path from the posterior
    [res,succ] = samplePosteriorNew(obs,A,states);
    if ~succ
       % It should be impossible to enter here
       %deb('Bad trace, resampling parameters');
       nResampled = nResampled + 1;
       continue;
    end
    trace = res.trace;
    
    % Store parameters and update posterior
    parSamples(n,:) = params;
    [aPost,bPost] = updateGamma(aPrior,bPrior, ...
            trace(:,2),trace(:,1), ...
            states,updates,rateFuncs);
    
    if (mod(n,PRINT_EVERY) == 0)
        fprintf('Taken %d samples.\n',n);
    end
    
    if mod(n,SAVE_EVERY) == 0
        filename = [outputFolder 'tempSamples' num2str(n)];
        fid = fopen(filename,'w');
        formatString = ['%f' repmat(' %f',1,R-1) '\n'];
        fprintf(fid,formatString,parSamples(1:n,:)');
        fclose(fid);
        fprintf('Wrote %d samples to file.\n',n);
    end
    
    n = n + 1;
end

%fprintf('Redrew parameters because of failure %d times.\n', nResampled);

end


function [new_a,new_b] = updateGamma(old_a,old_b,states,times,statespace,updates,rateFuncs)
% remember a, b are vectors
new_a = old_a;
new_b = old_b;

vecs = decodeStates(states,statespace);
rInd = findReactions(vecs,updates);
dt = diff(times);

nR = length(rateFuncs); % size(updates,1);

for r = 1:nR
    new_a(r) = new_a(r) + nnz(rInd == r);
    
    props = rateFuncs{r}(vecs(1:end-1,:));
    new_b(r) = new_b(r) + props' * dt;
end

end

function s = sampleGamma(a,b)
% Sample a number from a Gamma distribution with shape a and rate b
s = gamrnd(a, 1 ./ b);
end

function rInd = findReactions(states,updates)

jumps = diff(states);
[~,rInd] = myIsMember(jumps,updates);
end