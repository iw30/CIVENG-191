function prob=formLP(graph) %Fill the input graph=graph.edges
g=graph;
f=graph(:,3); %f Matrix
Aeq=[];
g1=graph(:,1);
g2=graph(:,2);
J=unique(g(:,1:2));
Z=zeros(1,length(g1));

% Building Aeq matrix
for i=1:(length(J)*2)
    if i<=length(J)
        Aeq(i,:)=Z;
        g1_1=find(g1==i);
        Aeq(i,g1_1)=1; 
    else
        Aeq(i,:)=Z;
        g2_1=find(g2==i-length(J));
        Aeq(i,g2_1)=1;
end
end

%Building beq matrix
beq=ones(size(Aeq(:,1)));

%Building lb and ub
lb=zeros(size(Aeq(1,:)))';
ub=ones(size(Aeq(1,:)))';

%Building Aineq and Bineq matrix

Aineq=[];
bineq=[];
if Aineq~=[]
bineq=zeros(size(A(1,:)))';
end

solver='linprog';
options=optimoptions('linprog','Algorithm','dual-simplex');
prob=struct('edges',g,'f',f,'Aineq',Aineq,'bineq',bineq,'Aeq',(Aeq),'beq',(beq),'lb',(lb),'ub',(ub),...
    'options',options,'solver',solver)
end


