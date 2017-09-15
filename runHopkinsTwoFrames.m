% Motion segmentation on John Hopkins155 dataset
% Pan Ji, pan.ji@adelaide.edu.au
% FEb 2017, University of Adelaide

clear, close all

addpath(genpath(pwd))
cd '/home/pan/Work/Hopkins155'; % cd to your Hopkins155 dir

file = dir;
idx = 0;
idx2 = 0;
idx3 = 0;
for i = 1:length(file)		
    if ( (file(i).isdir == 1) && ~strcmp(file(i).name,'.') && ~strcmp(file(i).name,'..') )
        filepath = file(i).name;
        eval(['cd ' filepath]);        
        
        f = dir;
        foundValidData = false;
        for j = 1:length(f)
            if ( ~isempty(strfind(f(j).name,'_truth.mat')) )
                ind = j;
                foundValidData = true;
				eval(['load ' f(ind).name]);
                break
            end
		end        
        cd ..        
        
        if(foundValidData && max(s)<=3)   
 			idx = idx+1;
 			n = max(s);           
            N = size(x,2);
            F = size(x,3);                        
            
            D = 3*F;
            X = reshape(permute(x(1:3,:,:),[1 3 2]),D,N); % data matrix with dimension DxF
            x1 = X(1:3,:);
            xlast = X(end-2:end,:);             
            
            X = zeros(9,N);
            for nn = 1:N
                X(:,nn) = kron(x1(:,nn),xlast(:,nn));
            end
            
            X = repmat(X,30,1);                       
			
            lambda1 = 0.23; lambda2 = 5.5; lambda3 = 1e5; kType = 'pol'; Affine = true; outlier = false; 
            param.alpha = 0.6;
            param.a = 2;
            param.b = 2; param.eta = 20;
            Missrate = lowRankKernelSubspaceClustering(X,s,lambda1,lambda2, lambda3, kType, Affine, outlier, param);
            
			disp([filepath ': ' num2str(100*Missrate) '%']);

 			missrate_tol(idx) = Missrate;
			if(max(s) == 2)
				idx2 = idx2+1;
				missrate_two(idx2) = Missrate;
			elseif(max(s) == 3)
				idx3 = idx3+1;
				missrate_three(idx3) = Missrate;
			end
 	
        end
    end
end
avgtol = mean(missrate_tol);
medtol = median(missrate_tol);
avgtwo = mean(missrate_two);
medtwo = median(missrate_two);
avgthree = mean(missrate_three);
medthree = median(missrate_three);
disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%.']);

