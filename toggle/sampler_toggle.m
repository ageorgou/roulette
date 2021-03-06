function parSamples = sampler_toggle(nSamples,obs,reactions,outputFolder)
%SAMPLER_TOGGLE Gibbs sampler with Russian Roulette truncation for the
%genetic toggle switch example.

global PRINT_EVERY SAVE_EVERY;

if nargin < 4 || isempty(outputFolder)
    outputFolder = './';
elseif outputFolder(end) ~= '/'
    outputFolder = [outputFolder '/'];
end

R = length(reactions); % number of reactions / parameters
parSamples = zeros(nSamples,R);

updates = cell2mat({reactions.update}');
rateFuncs = {reactions.propensity};
aPrior = [reactions.aPrior];
bPrior = [reactions.bPrior];

% Assume that initially gene1 is off and gene2 is on (we need to choose a
% concrete initial state so must somehow extend the initial observation)
initState = [obs(1,2:end) 0 1];
lowerBound = max(obs(:,2:end),[],1);

% Sample initial parameters...
oldPars = sampleGamma(aPrior,bPrior);

% ...initial truncation...
k = roulette();
limit = lowerBound + k;
oldSpace = makeStatespace(obs(:,2:end),updates,limit);
% ...and a trace using the initial truncation:
oldGen = makeGeneratorRest(oldSpace,updates,rateFuncs,oldPars);
succ = false;
while ~succ
    deb('Attemping to take first trace');
    [oldRes,succ] = samplePosteriorNew_toggle(obs,oldGen,oldSpace);
end
deb('Took first trace');
oldFull = oldRes.fullTrace;
oldTimes = oldFull(:,1);
oldTrace = oldFull(:,2:end);

aPost = aPrior; bPost = bPrior;
nAccept = 0;

for nS = 1:nSamples
    % Draw parameters from conditional posterior
    newPars = sampleGamma(aPost,bPost);
    
    % Choose a truncation level
    k = roulette();
    
    % Adjust the state-space
    lowerBound = max(obs(:,2:end),[],1);
    limit = lowerBound + k;
    newSpace = makeStatespace(obs(:,2:end),updates,limit);
    newGen = makeGeneratorRest(newSpace,updates,rateFuncs,newPars);
    
    oldGen = makeGeneratorRest(oldSpace,updates,rateFuncs,newPars);
    [oldGen,~] = uniformise(oldGen); % only need the discrete-time of this
    
    % Draw a trace (following Rao-Teh)
    [newRes,succ] = samplePosteriorNew_toggle(obs,newGen,newSpace);
    if ~succ
        deb('Bad trace, resampling parameters');
        continue;
    end
    deb('Taken new trace (%d)',nS);
    newFull = newRes.fullTrace;
    newTimes = newFull(:,1);
    newTrace = newFull(:,2:end);
        
    % Calculate the various likelihoods
    probNewNew = newRes.prob;
    
    [newGen,~] = uniformise(newGen); % now we need the discrete-time version
    oldDec = decodeStates(oldTrace,oldSpace);
    [oldReenc,found] = encodeStates(oldDec,newSpace);
    if all(found)
        succ = false;
        while ~succ
            deb('Calculating oldNew');
            [succ,probOldNew,~] = ... % oldS, new statespace
                FFBS_toggle(newGen,oldTimes,initState,obs,newSpace,oldReenc);
        end
    else
        probOldNew = 0;
        deb('probOldNew = 0');
    end
    
    newDec = decodeStates(newTrace,newSpace); % actual state vectors
    [newReenc,found] = encodeStates(newDec,oldSpace);
    if all(found)
        succ = false;
        while ~succ
            deb('Calculating newOld');
            [succ,probNewOld,~] = ... % newS, old statespace
                FFBS_toggle(oldGen,newTimes,initState,obs,oldSpace,newReenc);
        end
    else
        probNewOld = 0;
        deb('probNewOld = 0');
    end
    
    succ = false;
    while ~succ
        deb('Calculating oldOld');
        [succ,probOldOld,~] = ... 
            FFBS_toggle(oldGen,oldTimes,initState,obs,oldSpace,oldTrace);
    end
    
    deb('Calculated probabilities (%d)',nS);
    % Compute acceptance ratio
    if probNewOld == 0
        accRatio = 1;
    else
        accRatio = (probNewNew * probOldNew) / (probOldOld * probNewOld);
    end
    

    % Accept or reject sample
    if rand <= accRatio
        parSamples(nS,:) = newPars;
        oldPars = newPars;
        oldSpace = newSpace;
        oldTrace = newTrace;
        oldTimes = newTimes;
        % update posterior over parameters
        [aPost,bPost] = updateGamma(aPrior,bPrior, ...
            newRes.trace(:,2),newRes.trace(:,1), ...
            newSpace,updates,rateFuncs);
        nAccept = nAccept + 1;
    else
        parSamples(nS,:) = oldPars;
    end
    
    if mod(nS,PRINT_EVERY) == 0
        fprintf('Taken %d samples. Acceptance ratio: %f\n', ...
            nS,nAccept/nS);
    end
    
    if mod(nS,SAVE_EVERY) == 0
        filename = [outputFolder 'tempSamples' num2str(nS)];
        fid = fopen(filename,'w');
        formatString = ['%f' repmat(' %f',1,R-1) '\n'];
        fprintf(fid,formatString,parSamples(1:nS,:)');
        fclose(fid);
        fprintf('Wrote %d samples to file.\n',nS);
    end
    
end

end




function stopPoint = roulette()
global ROUL_CALLS ROUL_TERMS
ROUL_CALLS = ROUL_CALLS + 1;
REDUCTION_FACTOR = 0.95; % for example
shouldContinue = true;
stopPoint = 1;
%cumul = 1;
contProb = 1 / REDUCTION_FACTOR; %so guaranteed to continue the first time
while (shouldContinue)
    contProb = contProb * REDUCTION_FACTOR;
    %cumul = [cumul cumul(end)*contProb];
    if rand > contProb
        return;
    end
    stopPoint = stopPoint + 1;
    ROUL_TERMS = ROUL_TERMS + 1;
end
end


function states = makeStatespace(init,jumps,limit)

M = size(init,2); % number of observed species
MU = size(jumps,2) - M; % number of unobserved species
R = size(jumps,1); % number of reactions
% Expand the initial (observed) states to include the unobserved elements,
% i.e. promoters
init = unique(init,'rows');
nObs = size(init,1);
repObs = reshape(repmat(init,1,2^MU)',M,nObs*2^MU)';
unobs = repmat([0 0; 0 1; 1 0; 1 1],nObs,1);
limitFull = [limit 1 1]; % accounting for binary promoters

states = [repObs unobs];
added = states;
done = false;
while ~done
    newTargets = reshape(repmat(added,1,R)',(M+MU),size(added,1)*R)' ...
        + repmat(jumps,size(added,1),1);
    %added = newTargets(all(bsxfun(@le,newTargets,limit),2),:);
    % discard states outside the limits
    newTargets(any(bsxfun(@gt,newTargets,limitFull),2),:) = [];
    newTargets(any(newTargets < 0, 2),:) = [];
    % discard duplicates
    newTargets = unique(newTargets,'rows');
    % discards states that are already in statespace
    %alreadyIn = myIsMember(newTargets,states);
    alreadyIn = ismember(newTargets,states,'rows');
    %newTargets(alreadyIn,:) = [];
    added = newTargets(~alreadyIn,:);
    if isempty(added)
        done = true;
    end
    states = [states; added];
end
% add state IDs
states = [(1:size(states,1))' states];
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


function deb(msg,args)
% Print debugging message (uncomment when debugging)

% if nargin < 2
%     fprintf([msg '\n']);
% else
%     fprintf([msg '\n'],args);
% end
end
