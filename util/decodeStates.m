function stateVectors = decodeStates(stateIDs,lookup)
% DECODESTATES Transform a sequence of state IDs to their respective
% vectors.
%
% vecs = encodeStates(id,L) returns the vectors of the states referred to
% in id.
%
% L is an Nx(M+1) matrix containing the mapping between IDs and state
% vectors, where N is the number of the states and M is their dimension
% (number of components). id is an Rx1 column vector of state IDs, and vecs
% is an RxM matrix containing the corresponding state vectors.
%
% See also: encodeStates

[ind,loc] = ismember(stateIDs,lookup(:,1),'rows');
%[ind,loc] = myIsMember(stateIDs,lookup(:,1));

if ~all(ind)
    fprintf(2,'Could not recognise all state IDs.\n');
    %stateVectors = lookup(loc(ind),2:end);
end
stateVectors = lookup(loc,2:end);

end
