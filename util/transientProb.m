function result = transientProb(gen,initState,endState,time)

initProb = zeros(1,size(gen,1));
initProb(initState) = 1;

%endProb = initProb * expm(gen*time);
% Using EXPMV:
endProb = expmv(time,gen',initProb');

result = endProb(endState);

end
