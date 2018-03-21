function [ avg_Silhouette ] = avgSilhouette( data,cluster_assignment, K )

% This function computes the average Silhouette of a given clustering for a given data.
%variable 'cluster_assignment' is Nx1vector with an integer number between 1
%and K to denote which cluster the corresponding data point assigned to.
     N=size(cluster_assignment,1);
     counter=zeros(K,1);
   
     for n=1:N
     i=cluster_assignment(n);
     counter(i)=counter(i)+1;
     end
     
     PointToClusters=zeros(N,K);
     InCluster_a= zeros(N,1);
     nextBest_b=zeros(N,1);
     avg_Silhouette_point=zeros(N,1);
     avg_Silhouette_Cluster=zeros(K,1);
     
     
     
for dataPoint1=1:N
    nextPoint=dataPoint1+1;
    k1=cluster_assignment(dataPoint1);
    for dataPoint2= nextPoint:N
        k2=cluster_assignment(dataPoint2);
        sq_dis=sum((data(dataPoint1,:)-data(dataPoint2,:)).^2);
        
        PointToClusters(dataPoint1,k2)=PointToClusters(dataPoint1,k2)+sq_dis;
        PointToClusters(dataPoint2,k1)=PointToClusters(dataPoint2,k1)+sq_dis;
           
    end
    PointToClusters(dataPoint1,:)=PointToClusters(dataPoint1,:)./counter(:,1)';
    PointToClusters(dataPoint1,k1)=(PointToClusters(dataPoint1,k1)*counter(k1,1))/(counter(k1,1)-1);
    InCluster_a(dataPoint1)=PointToClusters(dataPoint1,k1);
    Sorted= sort(PointToClusters(dataPoint1,:));
    nextBest_b(dataPoint1)= Sorted(2);
    avg_Silhouette_point (dataPoint1)=(nextBest_b(dataPoint1)-InCluster_a(dataPoint1))/max(nextBest_b(dataPoint1),InCluster_a(dataPoint1));
    avg_Silhouette_Cluster(k1)=avg_Silhouette_Cluster(k1)+ avg_Silhouette_point (dataPoint1);
     
    
end

avg_Silhouette_Cluster=avg_Silhouette_Cluster./counter(:,1);

avg_Silhouette=sum(avg_Silhouette_Cluster)/K;

end

