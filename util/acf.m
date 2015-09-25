function rho = acf(x,maxLag)
%ACF Autocorrelation function for multi-dimensional series
%   ACF(X) where X is an NxM matrix, returns an NxM matrix R with the
%   autocorrelations of each column of X: R(i,j) is the autocorrelation of
%   X(:,j) at lag i-1.
%
%   ACF(X,L) returns an NxL matrix with the autocorrelations up to lag L.

if nargin < 2 || maxLag >= size(x,1)
    maxLag = size(x,1) - 1;
end

nCols = size(x,2);
rho = zeros(maxLag+1,nCols);

% Subtract mean - wrong?
y = bsxfun(@minus,x,mean(x));
%y = x;

for l = 0:maxLag
   % should technically be ... / size(x,1), but will cancel out with yVar
   rho(l+1,:) = bsxfun(@dot, y(1+l:end,:), y(1:end-l,:));
end

% Normalize
%yVar = y' * y;
%rho = rho / yVar;
yVar = bsxfun(@dot,y,y); % i.e. rho(1,:); should technically be / size(x,1)
rho = bsxfun(@rdivide,rho,yVar);

end