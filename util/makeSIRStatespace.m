function states = makeSIRStatespace(initState)

N = sum(initState);
N_S = initState(1);

nStates = (N+1)*(N+2)/2 - (N-N_S)*(N-N_S+1)/2;
states = zeros(nStates,length(initState) + 1);

n = 1;
while (N_S >= 0)
    N_rem = N - N_S;
    N_I = N_rem;
    while N_I >= 0
        N_R = N_rem - N_I;
        states(n,:) = [n N_S N_I N_R];
        N_I = N_I - 1;
        n = n + 1;
    end
    N_S = N_S - 1;
end

end