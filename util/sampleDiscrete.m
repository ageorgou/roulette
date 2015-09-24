function index = sampleDiscrete(probVector)

p = cumsum(probVector / sum(probVector));
u = rand;
index = find(p >= u, 1);

end


