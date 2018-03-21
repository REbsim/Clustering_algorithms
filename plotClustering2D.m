function [  ] = plotClustering2D(data,cluster_assignment,centerPoints)

% plotting 2-d dataset with colored clusters on
cols = {'r','g','b','y','c','k', 'm',[1,0.9,0.7],[0.8,0.7,0.7]};
K=size(centerPoints,1);
[ro co]=size(data);
 
if(co==2)
hold off;
    
for dataPoint = 1:ro
  xlabel('x1');
  ylabel('x2');   
  
    
   plot(data(dataPoint,1),data(dataPoint,2),'ko','markerfacecolor',cols{ cluster_assignment(dataPoint) });
   hold on
end



end




























end

