%solution should be a matrix of structures with elements indexed by year
%and condition, e.g.solution(year,condition+1). Each structure in the matrix should
%have a value field, and a decision field. For example,
%solution(3,3).value = 4, solution(3,3).decision = 'do nothing'.
%When you load transitionMatrices.mat it will give you an array of structures with fields matrix and
%decision. 
% In this lab, we assign transitionMatrices by
% load transitionMatrices
% cxu is a function to compute the cost function.
% Thus you might call the function below by typing solution = viteration ('transitionMatrices',@cxu,5) 
function solution = viteration(tmatrixFilename,cxuHandle,horizon)
% Your code starts here

load(tmatrixFilename) % Load tmatrixFilename
value=zeros(horizon,10); %Empty value matrix (matrix of zeros 5x10)
decision=repmat({''},horizon,10); %Empty decision matrix 5x10

for i=1:horizon % for loop for year 1 until 5 (backward)
    for j=0:9 %for loop condition 0-9 
        if i==1 % Assume V(5,j)=0
        % problem when the action is 'replacement'
        minp_1=cxu(j,'replacement'); 
        % problem when the action is 'maintenance'
        minp_2=cxu(j,'do nothing');
        % problem when the action is 'do nothing'
        minp_3=cxu(j,'maintenance');
        else  
        % problem when the action is 'replacement'
        minp_1=cxu(j,'replacement')+transitionMatrices(1).matrix(j+1,:)*value(end+2-i,:)'; 
        % problem when the action is 'maintenance'
        minp_2=cxu(j,'do nothing')+transitionMatrices(2).matrix(j+1,:)*value(end+2-i,:)' ; 
        % problem when the action is 'do nothing'
        minp_3=cxu(j,'maintenance')+transitionMatrices(3).matrix(j+1,:)*value(end+2-i,:)' ;
        end
        % Formulating Bellman equation and assigning the solution to value
        % matrix
        [value(end+1-i,j+1),idx]=min([minp_1,minp_2,minp_3]); %V(i,j)=min(...)
        if idx==1
            decision(end+1-i,j+1)={'replacement'};
        elseif idx==2
            decision(end+1-i,j+1)={'do nothing'};
        elseif idx==3
            decision(end+1-i,j+1)={'maintenance'};
        end
    end
end

%Final answer
solution.value=value;
solution.decision=decision;


% You code ends here
end