function [missrate, A, grps, obj, resid] = lowRankKernelSubspaceClustering(X, s, lambda1, lambda2, lambda3, kType, affine, outlier, param)
% Author: Pan Ji, University of Adelaide
% All rights reserved!
if(nargin<8)
    outlier = false;
end
if(nargin<7)
    affine = false;
end
if(nargin<6)
    kType = 'pol';
end
if(nargin<5)
    lambda3 = 1;
end
if(nargin<4)
    lambda2 = 1;
end
if(nargin<3)
    lambda1 = 1;
end

% normalize data to lie in [-1 1]
if(max(X(:)) <= 1 && min(X(:)) >= -1)
else
    X = X - min(X(:));
    X = X/(max(X(:))+1);
    X = 2*(X - 0.5); 
end

[~, N] = size(X);
alpha = param.alpha;

if(strcmpi(kType,'pol')) 
    a = param.a;
    b = param.b;
    KG = polynKernelMatrix(X,a,b); % hopkins 2.2 3   %two frames 3.4 2  % yale face 12 2    
elseif(strcmpi(kType,'lin'))
    KG = polynKernelMatrix(X,0,1);
elseif(strcmpi(kType,'rbf'))
    sig = std(sqrt(sum(X.^2)));
    KG = rbfKernelMatrix(X, 1*sig);    
else
    disp(['Select the kernel types: pol, lin, or rbf! Pol is used by default.']); 
    KG = polynKernelMatrix(X,1.4,3);
end


% ADMM parameters
epsilon = 1e-6; rho = 1e-8; maxIter = 1e3; eta = param.eta; max_rho = 1e10;

% Initializations
[U,S,V] = svd(KG);
B = U*sqrt(S)*V';
K = KG;

obj = [];
resid = [];
if(~outlier)
    if(~affine)        
        Y2 = zeros(N, N);
        A = zeros(N,N);
        iter= 0;
        
        while(iter<maxIter)
            iter = iter+1;
            
            % Update C
            tmp = A + Y2/rho;
            C = max(0, abs(tmp)-lambda1/rho) .* sign(tmp);
            C = C - diag(diag(C));
            
            % Update A
            K = B'*B;
            lhs = lambda2*K + rho*eye(N);
            rhs = lambda2*K - Y2 + rho*(C-diag(diag(C)));
            A = lhs\rhs;                       
            
            % Update B            
            tmp = KG-lambda2*(eye(N)-2*A'+A*A')/(2*lambda3);
            B = solveB(tmp, lambda3);            
           
            leq2 = A - (C - diag(diag(C)));
            stpC = max(abs(leq2(:)));
            if(iter == 1 || mod(iter,50)==0 || stpC<epsilon)
                disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ',stopALM=' num2str(stpC,'%2.3e')]);
            end
            
            if(stpC<epsilon)
                break;
            else                
                Y2 = Y2 + rho*leq2;
                rho = min(max_rho,rho*eta);
            end
            
        end
    else %affine        
        Y2 = zeros(N, N);
        y3 = zeros(1,N);
        A = zeros(N,N);
        iter= 0;
        
        while(iter<maxIter)
            iter = iter+1;
            
            % Update C
            tmp = A + Y2/rho;
            C = max(0, abs(tmp)-lambda1/rho) .* sign(tmp);
            C = C - diag(diag(C));
            
            % Update A
            K = B'*B;
            lhs = lambda2*K + rho*eye(N) + rho*ones(N,N);
            rhs = lambda2*K - Y2 - ones(N,1)*y3 + rho*(C-diag(diag(C))+ones(N,N));
            A = lhs\rhs;               
            
            % Update B            
            tmp = KG-lambda2*(eye(N)-2*A'+A*A')/(2*lambda3);
            B = solveB(tmp, lambda3);                     
            
            leq2 = A - (C - diag(diag(C)));
            leq3 = sum(A) - ones(1,N);
            
            obj(iter) = sum(svd(B)) + lambda1*sum(abs(C(:))) + ...
                0.5*lambda2*trace((eye(N)-2*A+A*A')*(B'*B)) + ...
                0.5*lambda3*norm(KG-B'*B,'fro')^2;
            resid(iter) = max(norm(leq2,'fro'),norm(leq3));
            
            stpC = max(abs(leq2(:)));
            stpC2 = max(abs(leq3));
            stpC = max(stpC, stpC2);
            if(iter == 1 || mod(iter,50)==0 || stpC<epsilon)                
                disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ',stopALM=' num2str(stpC,'%2.3e')]);
            end
            
            if(stpC<epsilon)
                break;
            else               
                Y2 = Y2 + rho*leq2;
                y3 = y3 + rho*leq3;
                rho = min(max_rho,rho*eta);                
            end
            
        end
    end
else% outliers
    if(~affine)
        %Y1 = zeros(N, N);
        Y2 = zeros(N, N);
        Y4 = zeros(N, N);
        A = zeros(N,N);
        E = zeros(N,N);
        iter= 0;
        
        while(iter<maxIter)
            iter = iter+1;
            
            % Update C
            tmp = A + Y2/rho;
            C = max(0, abs(tmp)-lambda1/rho) .* sign(tmp);
            C = C - diag(diag(C));
            
            % Update A
            K = B'*B;
            lhs = lambda2*K + rho*eye(N);
            rhs = lambda2*K - Y2 + rho*(C-diag(diag(C)));
            A = lhs\rhs;
            
            % Update B            
            tmp = KG - E - (0.5*lambda2*(eye(N)-2*A'+A*A')-Y4)/rho;
            B = solveB(tmp, rho);
            
            % Update E
            tmp = KG - B'*B + Y4/rho;
            E = max(0,tmp - lambda3/rho)+min(0,tmp + lambda3/rho);
            %E = solve_l1l2(tmp,lambda3/rho);                     
                       
            leq2 = A - (C - diag(diag(C)));
            leq4 = KG - B'*B - E;
            stpC =  max(abs(leq2(:)));
            stpC2 = max(abs(leq4(:)));
            stpC = max(stpC, stpC2);
            if(iter == 1 || mod(iter,50)==0 || stpC<epsilon)
                disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ',stopALM=' num2str(stpC,'%2.3e')]);
            end
            
            if(stpC<epsilon)
                break;
            else                
                Y2 = Y2 + rho*leq2;
                Y4 = Y4 + rho*leq4;
                rho = min(max_rho,rho*eta);
            end
            
        end
    else
        %Y1 = zeros(N, N);
        Y2 = zeros(N, N);
        y3 = zeros(1, N);
        Y4 = zeros(N, N);
        A = zeros(N, N);
        E = zeros(N, N);
        iter= 0;
        
        while(iter<maxIter)
            iter = iter+1;
            
            % Update C
            tmp = A + Y2/rho;
            C = max(0, abs(tmp)-lambda1/rho) .* sign(tmp);
            C = C - diag(diag(C));
            
            % Update A
            K = B'*B;
            lhs = lambda2*K + rho*eye(N) + rho*ones(N,N);
            rhs = lambda2*K - Y2 - ones(N,1)*y3 + rho*(C-diag(diag(C))+ones(N,N));
            A = lhs\rhs;            
            
            % Update B
            %tmp = 0.5*(K - E + KG + (Y1+Y4)/rho);
            tmp = KG - E - (0.5*lambda2*(eye(N)-2*A'+A*A')-Y4)/rho;
            B = solveB(tmp, rho);
            
            % Update E
            tmp = KG - B'*B + Y4/rho;
            E = max(0,tmp - lambda3/rho)+min(0,tmp + lambda3/rho);
            %E = solve_l1l2(tmp,lambda3/rho);
                        
            leq2 = A - (C - diag(diag(C)));
            leq3 = sum(A) - ones(1,N);
            leq4 = KG - B'*B - E;
            
            obj(iter) = sum(svd(B)) + lambda1*sum(abs(C(:))) + ...
                0.5*lambda2*trace((eye(N)-2*A+A*A')*(B'*B)) + ...
                lambda3*sum(abs(E(:)));
            resid(iter) = max(max(norm(leq2,'fro'),norm(leq3)), norm(leq4,'fro'));
            
            stpC = max(abs(leq2(:)));
            stpC2 = max(max(abs(leq3)), max(abs(leq4(:))));
            stpC = max(stpC, stpC2);
            if(iter == 1 || mod(iter,50)==0 || stpC<epsilon)
                disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ',stopALM=' num2str(stpC,'%2.3e')]);
            end
            
            if(stpC<epsilon)
                break;
            else
                %Y1 = Y1 + rho*leq1;
                Y2 = Y2 + rho*leq2;
                y3 = y3 + rho*leq3;
                Y4 = Y4 + rho*leq4;
                rho = min(max_rho,rho*eta);
            end
            
        end
    end
end

A = BuildAdjacency(thrC(C,alpha));
grp = SpectralClustering(A, max(s));
grps = bestMap(s,grp);
missrate = sum(s(:) ~= grps(:)) / length(s);

end

function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end
end































