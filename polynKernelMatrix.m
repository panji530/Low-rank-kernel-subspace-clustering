function K = polynKernelMatrix(X,c,d)
if(nargin<3)
	d = 2;
end
if(nargin<2)
	c = 0;
end

if(d == -2)
	tmp = X'*X;
	K   = double(tmp>0.9)+double(tmp<-0.9);
	%K   = K/trace(K);
elseif(d == -1)
	tmp = X'*X;
	K   = double(tmp>0.9)-double(tmp<-0.9);
	%K   = K/trace(K);
else
	K = (X'*X+c).^d;
	%K = K/trace(K);
end

end
