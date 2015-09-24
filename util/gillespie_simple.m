function trace = gillespie_simple(reactions,initState,stopTime)

rateFuncs = {reactions.rateFunc};

time = 0;
state = initState;
trace = [time state];

while (time < stopTime)
    % get rates from current state
    rates = cellfun(@(f) f(state),rateFuncs);
    
    % choose jump time and reaction
    R = sum(rates);
    % if no reactions are active, stop:
    if R == 0
        return
    end
    %rTime = exprnd(1/R);
    rTime = - log(rand) / R;
    x = cumsum(rates / R);
    u = rand;
    jumpID = find(x >= u, 1);
        
    % jump
    time = time + rTime;
    state = state + reactions(jumpID).update;
    trace = [trace; time state];
   
end