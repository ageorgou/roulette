function trace = gillespie_simple(reactions,initState,stopTime)
% GILLESPIE_SIMPLE Gillespie's Stochastic Simulation Algorithm.
%   GILLESPIE_SIMPLE(R,S,T) takes a structure of reactions R, an initial
%   state vector S and a final time T, and returns a sample trajectory from
%   the system defined in R. Each reaction must have two fields: "update",
%   a row vector containing the jump vector of the reaction; and
%   "rateFunc", a function handle representing its rate law as a function
%   of the state. S must be a row vector and its length N must match that
%   of the jump vectors. The output is a matrix with N+1 columns, where the
%   first column holds the jump times and the others hold the target state
%   vectors.

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