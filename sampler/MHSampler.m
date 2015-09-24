function parSamples = MHSampler(nSamples, obs, reactions)

global PRINT_EVERY

R = length(reactions); % number of reactions / parameters
parSamples = zeros(nSamples,R);
nAccept = 0;

updates = cell2mat({reactions.update}');
rateFuncs = {reactions.propensity};
aPrior = [reactions.aPrior];
bPrior = [reactions.bPrior];

%rouletteProbs = exp(-1:20);
%cumulProbs = [1 cumprod(rouletteProbs)];

maxStops = zeros(size(obs,1)-1,1); % what's the maximum number of terms we have encountered for each interval?
maxSpaces = {}; % the corresponding state space

oldPars = sampleGamma(aPrior,bPrior);
%oldPars = [2.5 0.01]; % for testing/debugging
oldPrior = gampdf(oldPars,aPrior,1./bPrior);
oldLL = loglike(obs,oldPars);
parSamples(1,:) = oldPars;

for nS = 2:nSamples
    % Propose parameters
    newPars = propose(oldPars);
    newPrior = gampdf(newPars,aPrior,1./bPrior);
    
    % Estimate the likelihood, working through pairs of observations
    newLL = loglike(obs,newPars);
    
    % Compute acceptance ratio (assume symmetric proposal)
    accRatio = exp(newLL - oldLL) * (newPrior / oldPrior);
    
    % Accept or reject sample
    if rand <= accRatio
        parSamples(nS,:) = newPars;
        oldPars = newPars;
        oldPrior = newPrior;
        oldLL = newLL;
        nAccept = nAccept + 1;
    else
        parSamples(nS,:) = oldPars;
    end
    
    if mod(nS,PRINT_EVERY) == 0
        fprintf('Taken %d samples. Acceptance ratio: %f\n', ...
            nS,nAccept/nS);
    end
end

fprintf('Acceptance ratio: %f\n',nAccept/nSamples);


    function LL = loglike(obs,pars)
        nObs = size(obs,1);
        times = obs(:,1);
        states = obs(:,2:end);
        LL = 0;
        %lowerBound = max(states,[],1); we want this for each segment
        
        for n = 2:nObs
            dt = times(n) - times(n-1);
            prob = 0;
            oldProb = 0;
            lowerBound = max(states(n-1,:), states(n,:));
            [stop,cumul] = roulette();
            initStates = states(n-1:n,:);
            for k = 1:stop
                % Check if we need to recompute the state-space
                if k > maxStops(n-1)
                    statespace = makeStatespace(initStates,updates,...
                        lowerBound+k);
                    maxStops(n-1) = k;
                    maxSpaces{n} = statespace;
                else
                   statespace = limitStatespace(maxSpaces{n},lowerBound+k);
                end
                [~,indices] = myIsMember(states(n-1:n,:),statespace(:,2:end));
                Q = makeGeneratorRest(statespace,updates,rateFuncs,pars);
                newProb = transientProb(Q,indices(1),indices(2),dt);
                % Prob. of remaining within the truncated state-space:
                %newProbAll = 1 - newSol(end);
                %newProb = newProbAll(encTarget) * newProbAll;
                prob = prob + (newProb - oldProb) / cumul(k+1);
                oldProb = newProb;
                initStates = statespace(:,2:end);
            end
            currLL = log(prob);
            LL = LL + currLL;
        end
        
    end

    function [stopPoint,cumul] = roulette()
        global ROUL_CALLS ROUL_TERMS
        ROUL_CALLS = ROUL_CALLS + 1;
        REDUCTION_FACTOR = 0.95; % for example
        shouldContinue = true;
        stopPoint = 1;
        cumul = 1;
        contProb = 1 / REDUCTION_FACTOR; %so guaranteed to continue the first time
        while (shouldContinue) 
            contProb = contProb * REDUCTION_FACTOR;
            cumul = [cumul cumul(end)*contProb];
            if rand > contProb
                return;
            end
            stopPoint = stopPoint + 1;
            ROUL_TERMS = ROUL_TERMS + 1;
        end
    end

end


function s = sampleGamma(a,b)
% Sample a number from a Gamma distribution with shape a and rate b
s = gamrnd(a, 1 ./ b);
end

function newSample = propose(oldSample)
global PROP_STD
newSample = oldSample + PROP_STD .* randn(size(oldSample));
while any(newSample < 0)
    newSample = oldSample + PROP_STD .* randn(size(oldSample));
end
end

function states = makeStatespace(init,jumps,limit)

M = size(jumps,2); % number of species
R = size(jumps,1); % number of reactions
states = init;
added = init;
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
    alreadyIn = ismember(newTargets,states,'rows');
    %alreadyIn = myIsMember(newTargets,states);
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
