function trace = gillespie(generator,initState,stopTime)

time = 0;
stateID = initState;
trace = [time stateID];

while (time < stopTime)
    % get rate from current state
    %rates = system.generator(stateID,:);
    rates = generator(stateID,:);
    rates(stateID) = 0; %ignore the diagonal entry
    
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
    stateID = find(x >= u, 1);
    
    % jump
    time = time + rTime;
    trace = [trace; time stateID];
    %trace = [trace; time system.states(stateID).vector];
   
end