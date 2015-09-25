function [I,loc] = myIsMember(A,B)
%MYISMEMBER Alternative to ISMEMBER(A,B,'rows')

dim = size(A,2);

% repeat each row of A for each row of B
nA = size(A,1);
nB = size(B,1);
repA = reshape(repmat(A,1,nB)',dim,nA*nB)';
repB = repmat(B,nA,1);

%all(bsxfun(@eq,repA,repB))
results = all(repA == repB,2);
results = reshape(results,nB,nA);
% now have one column per row of A

I = any(results,1)';
[indI, indJ] = find(results);
loc = zeros(nA,1);
loc(indJ) = indI;

end