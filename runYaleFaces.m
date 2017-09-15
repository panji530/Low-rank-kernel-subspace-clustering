% Face clustering on Yale face extended B dataset
% Pan Ji, pan.ji@adelaide.edu.au
clear, close all

load YaleBCrop025.mat Y;

nSet = [10 15 20 25 30 35 38];

for i = 1:length(nSet)	
    n = nSet(i);     
    for j = 1:(38-n+1)		
        X = [];
        s = [];
        for p =j:j+n-1
           X = [X Y(:,:,p)];
           s = [s,(p-j+1)*ones(1,64)];
        end        
		   
        [D,N] = size(X);     		
                
		lambda1 = 1.1e3; lambda2 = 2e-2; lambda3 = 1e5; kType = 'pol'; Affine = true; outlier = false; 
        param.alpha = 0.2;
        param.a = 12; 
        param.b = 2;
        param.eta = 20;
        [missrate, CMat, grp, obj, resid] = lowRankKernelSubspaceClustering(X,s,lambda1,lambda2, lambda3, kType, Affine, outlier,param);
						
        missrateTot{n}(j) = missrate;			
		
		disp([num2str(n) ' subjects, ' 'sequence ' num2str(j) ': ' num2str(100*missrateTot{n}(j)) '%']);		
	end	
    avgmissrate(n) = mean(missrateTot{n});
    medmissrate(n) = median(missrateTot{n});	
	disp([num2str(n) ' subjects: ']);
	disp(['Mean: ' num2str(100*avgmissrate(n)) '%, ' 'Median: ' num2str(100*medmissrate(n)) '%']);
end


