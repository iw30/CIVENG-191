% Dijkstra's algorithm for fidnings the shortest path from
% an origin node to all otehr nodes. The function takes an edgelist 
% of the form [node1 node2 edgeCost].

function [dist,prev] = myDijkstra(edges, origin)
%your code here

% dist[source]<- 0
dist=[];
dist(origin)=0;

% create vertex priority queue Q
pq=PriorityQueue;

graph2=edges(:,1:2);
uni=unique(reshape(graph2,[1,numel(graph2)]));
lgraph=numel(uni); %how many nodes in graph.edges

% for each vertex v in a graph where v is not source
dist(2:lgraph)=inf;
prev=zeros(1,lgraph);

% Q.add_with_priority(v,dist[v])

for i=1:lgraph
    pq.Insert(i,dist(i));
end

% while loop when Q is not empty

for z=1:lgraph
    u=pq.ExtractMin();
    lfind=find(edges(:,1)==u)'; %finding the first column of graph.edges=u
    edgesN=edges([lfind],:); %Simplifying graph.edges that contains u
    edgesN2=edgesN(:,1:2); %Simplified graph.edges two first columns
    edgesN3=edgesN(:,3); %Simplified graph.edges for the third column
    neighbor=reshape(edgesN2,[1,numel(edgesN2)]);
    neighbor(find(neighbor==u))=[]; 
    neighbor=sort(neighbor); %finding the neighboor of u in ascending orders
    
    for j=1:length(neighbor)
        lengthh=edgesN3(find(edgesN2(:,2)==neighbor(j)));
        alt=dist(u)+lengthh;
        if alt<dist(neighbor(j))
            dist(neighbor(j))=alt;
            prev(neighbor(j))=u;
            pq.DecreaseKey(neighbor(j),alt)
        end
    end
end
dist=dist;
prev=prev;
end