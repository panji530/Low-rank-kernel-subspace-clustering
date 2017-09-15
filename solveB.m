function B = solveB(K, rho)
% solve the following problem
% \min_B \|B\|_* + rho/2 \|B^TB - K\|_F^2
% Author: Pan Ji, University of Adelaide

K = 0.5*(K+K');
N = length(K);
[U, S, ~] = svd(K);
S = diag(S);

Lambda = zeros(N,1);
for i = 1:N
    si = S(i);
    l = roots([1 0 -si 1/(2*rho)]);
    l = l.*(imag(l)==0); % take real solutions
    l = max(l,0); % take positive solution
    l(length(l)+1) = 0;
    fval = 0.5*rho*(si*ones(size(l)) - l.^2).^2 + l;
    [~,idx] = min(fval);
    Lambda(i) = l(idx);
end

Lambda = diag(Lambda);
B = Lambda*U';

end
