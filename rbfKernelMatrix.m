function K = rbfKernelMatrix(X,sig)
if(nargin<2)
	sig = 1;
end
[~,N] = size(X);
K = X'*X/sig^2;
d = diag(K);
K = K - d*ones(1,N)/2 - ones(N,1)*d'/2;
K = exp(K);
end
