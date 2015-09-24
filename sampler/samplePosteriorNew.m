function [result,success] = samplePosteriorNew(obs,A,states)

% A is the infinitesimal generator matrix. states is a matrix each row of
% which corresponds to a state; the first column column contains state id's
% and the remaining columns the vector corresponding to that state.
% obs is a matrix of observations (time and measurement).


%updates = [reactions.update]';
%rateFuncs = {reactions.propensity};
% initState = ...

%states = makeStatespace(limits);
tF = obs(end,1);

initState = encodeStates(obs(1,2:end),states);
%for n = 2:nSamples
    
    % Sample a trace using the SSA
    sampleTrace = gillespie(A,initState,tF);
    %[S,T] = splitTrace(sampleTrace);
    S = sampleTrace(:,2:end);
    T = sampleTrace(:,1);
    % Uniformise the generator matrix
    [B,exitRate] = uniformise(A);
    
    % Sample from the subordinating process
    %U = virtualJumps(A,exitRate,S,T);
    % and merge the traces
    %[V,W] = merge(S,T,U);
    % Or all at once:
    [V,W] = addVirtualJumps(A,exitRate,S,T);
    
    % Run FFBS algorithm to get a new trace with times W
    success = false;
    while ~success
        deb('Running FFBS');
        %[succ,probV,newV] = FFBS(B,W,obs(1,2:end),obs(2:end,:),states,[]);
        [succ,probV,newV] = FFBS_noisy(B,W,obs(1,2:end),obs(2:end,:),states,[]);
        if ~succ
            success = false;
            result = [];
            deb('Need new parameters');
            return;
        end
        % Make sure the trace stays within the statespace and does not hit
        % the absorbing state (if any)
        if ~isWithinStatespace(newV,states)
            %success = false;
            deb('Outwith statespace');
        else
            success = true;
        end
    end
    
    % Drop the virtual jumps to obtain a new trace from the posterior process
    [newS,newT] = dropVirtual(newV,W);
    newTrace = [newT newS];
    
    
    result.trace = newTrace;
    result.fullTrace =[W newV];
    result.prob = probV;

end


function [newStates,newTimes] = dropVirtual(states,times)
jumpRows = [true; logical(diff(states))];
%alternatively:
% jumpRows = [1; find(diff(states)) + 1];
newStates = states(jumpRows);
newTimes = times(jumpRows);
end

function [newStates,newTimes] = addVirtualJumps(Q,exitRate,states,times)

newTimes = [];
newStates = [];
for i = 1:length(times)-1
   state = states(i);
   selfRate = Q(state,state) + exitRate;
   currTime = times(i);
   newTimes = [newTimes; currTime];
   newStates = [newStates; state];
   while currTime < times(i+1)
       % sample a new time from an exponential distribution
       sample = -log(rand) / selfRate;
       currTime = currTime + sample;
       if currTime < times(i+1)
           newTimes = [newTimes; currTime];
           newStates = [newStates; state];
       end
    end
end
newTimes = [newTimes; times(end)];
newStates = [newStates; states(end,:)];
end


function within = isWithinStatespace(stateIDs,statespace)
nStates = size(statespace,1);
within = ~any(stateIDs > nStates);
end

% function neighbours = neighbouringStates(stateVecs,centreState,maxDist)
% distances = abs(stateVecs - repmat(centreState,size(stateVecs,1),1));
% neighbours = sum(distances,2) <= maxDist;
% end

function deb(msg,args)
% Print debugging message (uncomment when debugging)

% if nargin < 2
%     fprintf([msg '\n']);
% else
%     fprintf([msg '\n'],args);
% end
end