function prob=branch(prob,graph) %fill prob=prob and graph=graph.edges
prob=prob;
graph=graph;

%Keep track what nodes are passed starting from node 1
    
    edges=prob.sol_edges;
    prob=[];
    Z=zeros(length(edges(:,1)),1);
    edges=[edges,Z];
    leftS=edges(:,1);
    allS=edges(find(leftS==1),:); % Start subtour progress from node 1
    rightS=allS(1,2);
    findS=find(edges(:,2)==rightS)
    edges(findS,3)=1 % Assign the zero equal to one if the edges is passed
    

% The while loop to keep track of the edges that are passed 
% until this reach node 1 again
    while rightS~=1      
        leftS=edges(:,1);
        allS=edges(find(leftS==rightS),:);
        rightS=allS(1,2);
        findS=find(edges(:,2)==rightS);
        edges(findS,3)=1;
    end

 % Now we can use the information from the third column matrix
  
 %Finding edges_solution with third column of 1
  edges11=edges(find(edges(:,3)==1),:); 
  edges12=edges11(:,1:2); % Delete the third column
  g1=graph(:,1:2); % Matrix graph.edges column 3 deleted
  member=ismember(g1,edges12,'rows'); % Comparing g1 with edges12
  g2=[graph,member]; %Adding the fourth column with member
  
fmember=find(member==1); %Finding what edge need to be branched

for i=1:length(fmember)
    myField=strcat('p',num2str(i));
    prob.(myField)=formLP(graph);
    Aineq=zeros(1,length(graph(:,1)));
    Aineq(fmember(i))=1;
    prob.(myField).Aineq=Aineq;
    prob.(myField).bineq=zeros(1,1);
    prob.(myField)=solveLP(prob.(myField));
end
end