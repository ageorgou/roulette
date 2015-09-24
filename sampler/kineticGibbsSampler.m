function paramSamples = kineticGibbsSampler(nSamples,obs,reactions,limits)

global PRINT_EVERY;

updates = cell2mat({reactions.update}');
rateFuncs = {reactions.propensity};
a = [reactions.aPrior];
b = [reactions.bPrior];
nPar = length(reactions); % number of parameters

% Make state-space
states = makeStatespace(obs(:,2:end),updates,limits);

paramSamples = zeros(nSamples,nPar);

% Draw parameters from the prior
currentParam = sampleGamma(a,b);
paramSamples(1,:) = currentParam;

n = 2;
nResampled = 0;
while n <= nSamples
    
    % Create the generator
    A = makeGeneratorRest(states,updates,rateFuncs,currentParam);
    
    % Draw a new path from the posterior
    [res,succ] = samplePosteriorNew(obs,A,states);
    if ~succ
        %deb('Bad trace, resampling parameters');
        nResampled = nResampled + 1;
        continue;
    end
    
    trace = res.trace;
    newT = trace(:,1);
    newS = trace(:,2:end);
    
    % Get the posterior over the parameters and sample from it
    [a_new,b_new] = updateGamma(a,b,newS,newT,states,updates,rateFuncs);
    currentParam = sampleGamma(a_new,b_new);
    paramSamples(n,:) = currentParam;
    if (mod(n,PRINT_EVERY) == 0)
        fprintf('Taken %d samples.\n',n);
    end
    n = n + 1;
end

fprintf('Redrew parameters because of failure %d times.\n', nResampled);

end

function states = makeStatespace(init,jumps,limit)

M = size(jumps,2); % number of species
R = size(jumps,1); % number of reactions
states = unique(init,'rows');
added = states;
done = false;
while ~done
    newTargets = reshape(repmat(added,1,R)',M,size(added,1)*R)' ...
        + repmat(jumps,size(added,1),1);
    %added = newTargets(all(bsxfun(@le,newTargets,limit),2),:);
    % discard states outside the limits
    newTargets(any(bsxfun(@gt,newTargets,limit),2),:) = [];
    newTargets(any(newTargets < 0, 2),:) = [];
    % discard duplicates
    newTargets = unique(newTargets,'rows');
    % discards states that are already in statespace
    alreadyIn = myIsMember(newTargets,states);
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