function [ bestClustering,bestMeans, bestCovs, bestPriors ] = GMM(data,K )
% this function implements Gaussian mixture model. It runs the model for
% 100 times and returns the best results acheived in terms of the highest
% lower bound B of the likelihood.
% Input: dataset data and number of clusters K
% Output: 
%variable 'bestClustering' is Nx1vector with an integer number between 1
%and K to denote which cluster the corresponding data point assigned to.
%variables 'bestMeans' and 'bestCovs' return the parametes of the best
%Gaussian built on the dataset in 100 runs.

[N f]=size (data);
bestClustering= zeros(N,1); 
bestMeans=zeros(K ,f);
bestCovs=zeros(f,f,K);
maxB=-inf;

for run=1:100
covs=zeros(f,f,K);
new_covs=zeros(f,f,K);
clustering=zeros(N,1);
priors =zeros(1,K);
MaxIts = 150;

B = [];
B(1) = -inf;
converged = 0;
q = zeros(N,K);
it = 1;
tol = 1e-2;
% initialization priors and gaussians

for i=1:f
    % helping GMM starting with good random means
a = min(data(:,i));
b = max(data(:,i));
means(1:K,i) = a + (b-a).*rand(K,1) ;
end
for k = 1:K
    priors(1,k)=1/K;
    covs(:,:,k) = rand*eye(f);
end


 while 1
    it = it + 1;
    % logp holds log (p(n,k))
    
    for k = 1:K
        const = -(f/2)*log(2*pi) - 0.5*log(det(covs(:,:,k)));
        for point=1:N
                 dif=data(point,:)-means(k,:);
                logp(point,k) = const - 0.5 * (dif*inv(covs(:,:,k))*dif');
         end
     
    end
              
    % Compute the current lower bound on the likelihood B
       if it>1
                
        term1=0;
        term2=0;
        term3=0;
        for k=1:K
           for n=1:N 
            term1=term1+q(n,k)*log(priors(k))';
            term2=term2+q(n,k)*logp(n,k)';
            term3=term3+q(n,k)*log(q(n,k)');
           end  
        end
              B(it) = term1+term2-term3;
                
                if abs(B(it)-B(it-1))<tol
                    converged = 1;

                end
            end

            if converged == 1 || it> MaxIts
                break
            end
   %updating q(n,k)
    logPriors=log( priors);
    logPPriors = logp + repmat(logPriors,N,1);
    qNominator = exp(logPPriors);
    % avoiding numerical problems    
    qNominator(qNominator<1e-3) = 1e-3;
    qNominator(qNominator>1-1e-3) =1-1e-3;
    
    q = qNominator./repmat(sum(qNominator,2),1,K);
    sumQ=sum(q);
        

% updating priors   
 priors= sumQ./N   ;

 % updatg means 
 for k = 1:K
     means(k,:)=0;
     for n = 1:N
         means(k,:) = means(k,:)+ q(n,k)*data(n,:);
     end
   means(k,:)  =means(k,:)./ sumQ(k);
 end
   
% updating covariances
 
 
 for k = 1:K
    new_covs(:,:,k)=zeros(f,f); 
 for n = 1:N
   Xm=data(n,:)-means(k,:);
 new_covs(:,:,k)= new_covs(:,:,k)+ (q(n,k).*Xm)'*(Xm);
 end
 covs(:,:,k)= new_covs(:,:,k)./sumQ(k);
 end
 
end  
for i=1:N
        [value index]= max(q(i,:));
       clustering(i)= index;
        
end
% Here, I get rid of the runs with number of clusters less than K
% for which Silhouete equals NAN.

  Sil= avgSilhouette(data,clustering, K);
 if ~isnan(Sil)
      if B(it)>maxB
          maxB=B(it);
          bestClustering= clustering ;
          bestMeans=means;
          bestCovs=covs;
          bestPriors=priors;
          fprintf('OK');
      end
 end
 
end

   end