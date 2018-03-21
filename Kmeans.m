function [ cluster_assignment cluster_centoids d ] = Kmeans( data,K , DistMeasure )
% this function implements K-means algorithm to partition a given dataset
%into k clusters according either with Manhatan 'Manh' or eucldian distance
%'Ecu'
%variable 'cluster_assignment' is Nx1vector with an integer number between 1
%and K to denote which cluster the corresponding data point assigned to.
%variables 'cluster_centoids' the centroids.

[N  f]= size(data);
IsConverged=false;
%cluster_centoids=rand(K,f)*10-4;
for i=1:f
    % helping Kmeans starting with good random means
a = min(data(:,i));
b = max(data(:,i));
cluster_centoids(1:K,i) = a + (b-a).*rand(K,1) ;
end





d=0;
% Initial assignment of datapoints to clustors. 
%cluster_assignment(i) holds the clustor number to which datapoint i is assigned.
cluster_assignment=zeros(N,1);
old_cluster_assignment=zeros(N,1);

% count(i) holds the number of datapoints in clustor i.
count=zeros(K,1);
 
while ~IsConverged
    
count=zeros(K,1);
%cluster assignment step
distance=zeros(N,K);
for dataPoint=1:N
    
    for cluster=1:K
    if (strcmpi (DistMeasure, 'Manh'))
    distance(dataPoint,cluster)=sum(abs(data(dataPoint,:)-cluster_centoids(cluster,:)));
        
    else
    distance(dataPoint,cluster)=sum((data(dataPoint,:)-cluster_centoids(cluster,:)).^2);
    end
    
    end
   
    
    [minDistance  closest_cluster]= min(distance(dataPoint,:));

    cluster_assignment(dataPoint)=closest_cluster;
    count(closest_cluster)= count(closest_cluster)+1;
    %d=d +minDistance;
    
end


% calculating new centroids step
new_centroids=zeros(K,f);

for dataPoint=1:N
        new_centroids(cluster_assignment(dataPoint),:)= new_centroids(cluster_assignment(dataPoint),:)+ data(dataPoint,:);
end

for i=1:K
 cluster_centoids(i,:)=(1/count(i)).*new_centroids(i,:);
end
 if (cluster_assignment==old_cluster_assignment)
     IsConverged= true;
 else 
     old_cluster_assignment=cluster_assignment;
 end
 
end




for dataPoint=1:N
       d= d+sqrt(sum((data(dataPoint,:)-cluster_centoids(cluster_assignment(dataPoint),:)).^2));
end





end

