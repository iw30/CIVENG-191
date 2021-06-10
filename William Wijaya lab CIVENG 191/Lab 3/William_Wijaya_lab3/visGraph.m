function s=visGraph(edges)
e_12=edges(:,1:2);
e3=edges(:,3);
C=unique(edges(:,1:2));
Z=zeros(1,length(C));
Com=[];
for i=1:length(C)
    Com(i,:)=Z;
    e1=edges(:,1);
    e1_1=e_12(find(e1==i),:);
    e2=e1_1(:,2);
    Com(i,e2)=e3(find(e1==i),:);
end
bg=biograph(Com,[],'ShowWeights','on');
view(bg)