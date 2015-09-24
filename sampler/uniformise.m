function [Q_disc,exitRate] = uniformise(Q_cont,testExitRate)
qMax = -min(diag(Q_cont));
if (nargin < 2 || testExitRate < qMax)
    exitRate = 1.1 * qMax;
else
    exitRate = testExitRate;
end

%Q_disc = Q_cont / exitRate + eye(size(Q_cont));
Q_disc = Q_cont / exitRate + speye(size(Q_cont));
end