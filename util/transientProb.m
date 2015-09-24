function result = transientProb(gen,initState,endState,time)

initProb = zeros(1,size(gen,1));
initProb(initState) = 1;

%endProb = initProb * expm(gen*time);
% Uniformisation testing:
%endProb2 = FoxGlynnTransient(gen,initProb',time,0.00001);
%endProb = FoxGlynnTransient(gen,initProb',time,0.00001);

% Using EXPMV:
endProb = expmv(time,gen',initProb');

result = endProb(endState);

end