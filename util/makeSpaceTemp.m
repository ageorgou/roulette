function states = makeSpaceTemp(init,jumps,limit)

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