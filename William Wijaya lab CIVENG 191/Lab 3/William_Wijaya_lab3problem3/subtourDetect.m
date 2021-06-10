function hasSubtours=subtourDetect(sol_edges)
edges=sol_edges;
Z=zeros(length(edges(:,1)),1);
edges=[edges,Z];
leftS=edges(:,1);
allS=edges(find(leftS==1),:); % Start subtour progress from node 1
rightS=allS(1,2);
findS=find(edges(:,2)==rightS);
edges(findS,3)=1; % Assign the zero equal to one if the edges is passed

% The while loop to keep track of the edges that are passed 
% until this reach node 1 again
while rightS~=1
    leftS=edges(:,1);
    allS=edges(find(leftS==rightS),:);
    rightS=allS(1,2);
    findS=find(edges(:,2)==rightS);
    edges(findS,3)=1;
end

% All of the above program will generate 3 column matrix.
% The first two column matrix will be the sol_edges
% The third column matrix is a matrix of 0 or 1

% If all third column matrix is equal to one, 
% then all edges starting from one are passed, meaning there is no
% subtours.
if all(edges(:,3))==1
    hasSubtours=0;
else
    hasSubtours=1;
end

