function H = ess(x)
%ESS Compute effective sample size.
%  ESS(X), where X is an NxM matrix, returns a 1xM vector containing the
%  effective sample size of X, computed separately across each of the M
%  dimensions.

N = size(x,1);
dim = size(x,2);

autocorr = acf(x); % this includes autocorrelation at lag 0, i.e. 1

%H = N ./ (1 + 2*sum(autocorr(2:end,:)));
% but we should stop once the autocorrelations become negative

H = zeros(1,dim);
for p = 1:dim
    if all(autocorr(:,p) > 0)
        stop = N;
    else
        stop = find(autocorr(:,p) <= 0,1);
    end
    H(p) = N / (1 + 2*sum(autocorr(2:stop,p)));
end

end