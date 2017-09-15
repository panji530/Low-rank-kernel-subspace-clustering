function K = polynKernelMatrix(X,c,d)
% Author: Pan Ji, University of Adelaide
if(nargin<3)
	d = 2;
end
if(nargin<2)
	c = 0;
end

if(d == -2)
	tmp = X'*X;
	K   = double(tmp>0.9)+double(tmp<-0.9);	
elseif(d == -1)
	tmp = X'*X;
	K   = double(tmp>0.9)-double(tmp<-0.9);	
else
	K = (X'*X+c).^d;	
end

end
