function tour=tourFun(sol_edges) %fill sol_edges=prob.sol_edges
trackTour=[];
edges=sol_edges;
leftS=edges(:,1);
allS=edges(find(leftS==1),:); % Start subtour progress from node 1
tracktour=[1];
rightS=allS(1,2);
tracktour(1,2)=rightS;

i=3;
while rightS~=1
    leftS=edges(:,1);
    allS=edges(find(leftS==rightS),:);
    rightS=allS(1,2);
    tracktour(1,i)=rightS;
    i=i+1;
end

tour=tracktour;
