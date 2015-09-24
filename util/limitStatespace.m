function newSpace = limitStatespace(space,limit)


toKeep = all(bsxfun(@le,space(:,2:end),limit),2);
%newSpace = [(1:nnz(toKeep))' space(toKeep,2:end)];
newStates = space(toKeep,2:end);
newSpace = [(1:size(newStates,1))' newStates];

end