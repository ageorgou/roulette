function [stateIDs,found] = encodeStates(stateVectors,lookup)
% ENCODESTATES Transform a sequence of state vectors to their respective
% IDs.
%
% IDs = encodeStates(vectors,L) returns the IDs of the states contained in
% vectors.
%
% L is an Nx(M+1) matrix containing the mapping between IDs and state
% vectors, where N is the number of the states and M is their dimension
% (number of components). Each row in the RxM matrix vectors represents one 
% state vector. The return value IDs is an Rx1 column vector whose elements
% are the IDs corresponding to the states in vectors.
%
% See also: decodeStates

%[ind,loc] = ismember(stateVectors,lookup(:,2:end),'rows');
[ind,loc] = myIsMember(stateVectors,lookup(:,2:end));

% if ~all(ind)
%     fprintf(2,'Could not find all states in statespace.');
% end
stateIDs = lookup(loc(ind),1);
found = ind;
% for absorbing states, represented as rows of NaN's (maybe not needed):
%NaN_rows = all(isnan(stateVectors),2);
%stateIDs = [stateIDs; lookup(NaN_rows,1)];

end
