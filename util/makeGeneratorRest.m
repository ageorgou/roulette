function Q = makeGeneratorRest(states,updates,rates,pars)
% MAKEGENERATORREST Create a sparse generator matrix for a given system,
% including an absorbing state.
%
% makeGenerator(states,updates,rates) returns a sparse matrix Q that
% represents the infinitesimal generator matrix for a system. states is the
% state-space of the system, represented as an array. updates is a matrix
% of the jump vectors, with one row for each reaction. rates is a cell
% array of function handles, giving the rate of each reaction as a function
% of the state. The resulting matrix Q is of size (N+1) x (N+1), where N is
% the size of the state-space (number of rows of the matrix states). The
% additional row andcolumn correspond to an absorbing state, representing
% states of the system outside the specified state-space.
%
% See also: makeStatespace

nStates = size(states,1);
nReacs = size(updates,1);
absState = nStates + 1;

startInd = [];
endInd = [];
jumpRates = [];
exitRates = zeros(nStates,1);
rest = zeros(nStates,1);

for n = 1:nReacs
    endStates = bsxfun(@plus,states(:,2:end),updates(n,:));
    %[found,newEndInd] = myIsMember(endStates,states(:,2:end));
    [found,newEndInd] = ismember(endStates,states(:,2:end),'rows');
    newStartInd = states(found,1);
    % for vectorised rate functions:
    newJumpRates = pars(n) * rates{n}(states(:,2:end));
    % otherwise:
%   newJumpRates = zeros(length(newStartInd),1);
%   for k = 1:length(newStartInd)
%     newJumpRates(k) = rates{n}(states(newStartInd(k),2:end));
%   end;
    startInd = [startInd; newStartInd];
    endInd = [endInd; newEndInd(found)];
    jumpRates = [jumpRates; newJumpRates(found)];
    rest(~found) = rest(~found) + newJumpRates(~found);
    exitRates = exitRates + newJumpRates;
end

startToRest = find(rest ~= 0);
restTarget = repmat(absState,length(startToRest),1);

Q = sparse([startInd; startToRest; (1:nStates)'], ...
    [endInd; restTarget; (1:nStates)'], ...
    [jumpRates; rest(startToRest); -exitRates], ...
    nStates+1,nStates+1);

end