function [cluster_assignment medoids distance] = Kmedoids( data,K )
%This function implements K-medoids algorithm by Park & Jun 2009 to partition a given dataset
%into k clusters. 
%variable 'cluster_assignment' is Nx1vector with an integer number between 1
%and K to denote which cluster the corresponding data point assigned to.

[n f]= size(data);

sequared_distance=zeros(n);
cluster_assignment=zeros(n,1);
medoids=zeros(K,f);
medoids_index=zeros(K,1);
V=zeros(n,2);
stop=false;

% calculating inter-point distances
for dataPoint1=1:n
    
    for dataPoint2=dataPoint1:n
        
    sequared_distance(dataPoint1,dataPoint2)=sqrt(sum((data(dataPoint1,:)-data(dataPoint2,:)).^2));
    sequared_distance(dataPoint2,dataPoint1)=sequared_distance(dataPoint1,dataPoint2);
    end
       
end
% medoids initialization
for Point_j=1:n
    
    V(Point_j,1)=0;
    V(Point_j,2)=Point_j;
    
    for Point_i=1:n
    denominator= sum(sequared_distance(Point_i,:));
    V(Point_j,1)= V(Point_j,1)+ sequared_distance(Point_i,Point_j)/denominator;
    end
end
SortedV= sortrows(V);
for k=1:K
  medoids(k,:)= data(SortedV(k,2),:);
  medoids_index(k,1)= SortedV(k,2); 
end  
%initialMedoids =medoids

% cluster Assignment
distance=0;
point_assignment=zeros(K,n);
cluster_counter=zeros(K,1);


for dataPoint1=1:n
      cluster_assignment(dataPoint1)=1;
      ShortestdistanceToMedoid=sequared_distance(dataPoint1,medoids_index(1,1));
  for k=2:K
      if (sequared_distance(dataPoint1,medoids_index(k,1))<ShortestdistanceToMedoid)
      cluster_assignment(dataPoint1)=k;
      ShortestdistanceToMedoid=sequared_distance(dataPoint1,medoids_index(k,1));
      end   
      
  end
   
distance=distance + ShortestdistanceToMedoid;
cluster=cluster_assignment(dataPoint1);
cluster_counter(cluster,1)= cluster_counter(cluster,1)+1;
point_assignment(cluster,cluster_counter(cluster,1))=dataPoint1;

end

while ~stop
% updating medoids 
% for each cluster find the point that has the total lowest distance from all other points in the cluster
intraDistance= zeros(K,n);

for k=1:K
    counter=cluster_counter(k);
  for point1=1:counter
    
      for point2=point1:counter
      temp=sequared_distance (point_assignment(k,point1), point_assignment(k,point2));
      intraDistance(k,point1) =intraDistance(k,point1)+temp;
      intraDistance(k,point2) =intraDistance(k,point2)+temp;    
      end
    
  end
   [value index]=min(intraDistance(k,1:counter)) ;
      
      
   medoids_index(k,1)= point_assignment(k,index);
   medoids(k,:) =data(medoids_index(k,1));
   
end

% cluster Assignment
new_distance=0;
point_assignment=zeros(K,n);
cluster_counter=zeros(K,1);


for dataPoint1=1:n
      cluster_assignment(dataPoint1)=1;
      ShortestdistanceToMedoid=sequared_distance(dataPoint1,medoids_index(1,1));
  for k=2:K
      if (sequared_distance(dataPoint1,medoids_index(k,1))<ShortestdistanceToMedoid)
      cluster_assignment(dataPoint1)=k;
      ShortestdistanceToMedoid=sequared_distance(dataPoint1,medoids_index(k,1));
      end   
      
  end

new_distance=new_distance + ShortestdistanceToMedoid;
cluster=cluster_assignment(dataPoint1);
cluster_counter(cluster,1)= cluster_counter(cluster,1)+1;
point_assignment(cluster,cluster_counter(cluster,1))=dataPoint1;

end

if (new_distance== distance)
    stop= true;
else 
 distance= new_distance;
end 


end



end

