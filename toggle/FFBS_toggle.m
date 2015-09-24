function [success,totalProb,sample] = FFBS_toggle(P,times,initState,obs,statespace,oldSample)

dim = size(P,1); % number of states
nTimes = length(times);
sample = zeros(size(times));
totalProb = 1;

initState = encodeStates(initState,statespace);

probs = zeros(dim,nTimes); % observation probabilities

% Calculate the forward messages (from the past)
a = zeros(dim,nTimes);
% Initialize by assuming we know (noiselessly) the first state
a(:,1) = zeros(dim,1); a(initState,1) = 1;
probs(initState,1) = 1;

for ii = 2:nTimes
    % Find which (if any) observations were made since the last time
    ind = obs(:,1) >= times(ii-1) & obs(:,1) < times(ii);
    if ~any(ind)
        probs(:,ii) = ones(dim,1);
    elseif nnz(ind) == 1
        probs(:,ii) = obsProbs(obs(ind,2:end),statespace);
    else
        %fprintf(2,'More than one observation per interval found, proceed with caution.\n');
        nObs = nnz(ind);
        probsInt = zeros(dim,nObs);
        O = obs(ind,2:end);
        for jj = 1:nObs
            probsInt(:,jj) = obsProbs(O(jj,:),statespace);
        end
        probs(:,ii) = prod(probsInt,2);
    end
    % Recursion for the forward message:
    %a(:,ii) = sum(repmat(a(:,ii-1) .* probs(:,ii-1),1,dim) .* P)';
    a(:,ii) = ((a(:,ii-1) .* probs(:,ii-1))' * P)';
    if ~any(a(:,ii))
        success = false;
        fprintf('Failed at forward step (%d).\n',ii);
        return;
    end
end
% Finish by computing the observation probabilities for the last time point
ind = obs(:,1) == times(nTimes);
if any(ind)
    probs(:,nTimes) = obsProbs(obs(ind,2:end),statespace);
end


% Calculate the backward messages, sampling a state at each time point
b = zeros(dim,nTimes);
b(:,nTimes) = a(:,nTimes) .* probs(:,nTimes);
if ~any(b(:,nTimes))
        success = false;
        return;
end

% Are we sampling a new trajectory or just calculating the probability of a
% given one (oldSample) ? In the latter case, we just copy the old
% trajectory and use its states in the probability calculation.
samplingNew = nargin < 6 || isempty(oldSample);
if samplingNew
    % Sample a final state according to b(:,nTimes)...
    randInd = sampleDiscrete(b(:,nTimes));
else
    randInd = oldSample(nTimes);
end
sample(nTimes) = randInd;
totalProb = totalProb * (b(randInd,nTimes) / sum(b(:,nTimes)));
% ...and repeat going back
for ii = nTimes-1:-1:1
    % Recursion for the backward message:
    b(:,ii) = a(:,ii) .* probs(:,ii) .* P(:,sample(ii+1));
    if samplingNew
        randInd = sampleDiscrete(b(:,ii));
    else
        randInd = oldSample(ii);
    end
    sample(ii) = randInd;
    totalProb = totalProb * (b(randInd,ii) / sum(b(:,ii)));
end
success = true;
%     Probably not worth defining / using in the above:
%     function [newState,prob] = handle_b(messageVec)
%        randInd = sampleDiscrete(messageVec);
%        newState = randInd;
%        prob = messageVec(randInd) / sum(messageVec);
%     end

end

function p = obsProbs(obsVector,statespace)
% Modified for the toggle switch case, i.e. where only some species are
% observed. obsVector is a 2x1 vector with observed values of the
% repressors at a particular time point.

D = 10^-6;

% Use the Euclidean distance from the observed state... (only taking into
% account the first two dimensions, i.e. the repressors, and disregarding
% the promoters)
diffs = bsxfun(@minus,statespace(:,2:3),obsVector);
dists = sqrt(sum(diffs.^2,2));
p = 1 ./ (2 .^ dists + D); % ...and geometric-like observation probability
p = [p; 0]; % Add an entry for the absorbing state (prob. 0 to be observed)
p = p / sum(p); % Normalize (maybe not actually required?)
end


% Cut off:

% totalProb = 1; % probability of the full trajectory sampled
% if isempty(oldSample)
%     % Sample a final state according to b(:,nTimes)...
%     randInd = sampleDiscrete(b(:,nTimes));
%     sample(nTimes) = randInd;
%     totalProb = totalProb * (b(randInd,nTimes) / sum(b(:,nTimes)));
%     % ...and repeat going back
%     for ii = nTimes-1:-1:1
%         % Recursion for the backward message:
%         b(:,ii) = a(:,ii) .* probs(:,ii) .* P(:,sample(ii+1));
%         randInd = sampleDiscrete(b(:,ii));
%         sample(ii) = randInd;
%         totalProb = totalProb * (b(randInd,ii) / sum(b(:,ii)));
%     end
%     success = true;
% else
%     % Calculate the probability of the last observation using b(:,nTimes)
%     totalProb = totalProb * ...
%         (b(oldSample(nTimes),nTimes) / sum(b(:,nTimes)));
%     % ...and repeat going back
%     for ii = nTimes-1:-1:1
%         % Recursion for the backward message:
%         b(:,ii) = a(:,ii) .* probs(:,ii) .* P(:,sample(ii+1));
%         randInd = sampleDiscrete(b(:,ii));
%         sample(ii) = randInd;
%         totalProb = totalProb * (b(randInd,ii) / sum(b(:,ii)));
%     end
%     success = true;
% end


