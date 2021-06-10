%solution should be an array of structures with elements indexed by year,
%e.g., solution(state). Each structure in the matrix should
%have a value field, and a decision field. For example,
%solution(3).value = 4, solution(3).decision = 'do nothing'.
%When you load transitionMatrices.mat it will give you an array of structures with fields matrix and
%decision. 
% In this lab, we assign transitionMatrices by
% 'load transitionMatrices'
% cxu is a function to compute the cost function.
% Thus you might call the function below by typing solution = viterationInf ('transitionMatrices',@cxu,0.001) 

function solution = viterationInf (tmatrixFilename, cxuHandle, epsilon)
% In the case of infinite horizon, solution is not a function over time,
% it is only a function of state, i.e.
% solution(state).value, solution(state).decision

%You code starts here

%Fill epsilon=0.001
load(tmatrixFilename) %Load tmatrixFile name
value=zeros(1,10); %Empty value matrix (matrix of zeros nx10, it acts like a history.
decision=repmat({''},1,10); %Empty decision matrix nx10
difference=inf; %Set difference=inf as a starting point
diffmat=repmat(difference,3,10); %Difference matrix
epsmat=repmat(epsilon,1,10); %Epsilon matrix
alpha=0.95; %Discount factor

i=1; 
while all(abs(diffmat(i+1,:)-diffmat(i,:))<epsmat(1,:))==0
    for j=0:9 %condition 0 until 9
        if i==1
        minp_1=cxu(j,'replacement'); 
        % problem when the action is 'maintenance'
        minp_2=cxu(j,'do nothing');
        % problem when the action is 'do nothing'
        minp_3=cxu(j,'maintenance');
        else  
        % problem when the action is 'replacement'
        minp_1=cxu(j,'replacement')+alpha*transitionMatrices(1).matrix(j+1,:)*value(i-1,:)'; 
        % problem when the action is 'maintenance'
        minp_2=cxu(j,'do nothing')+alpha*transitionMatrices(2).matrix(j+1,:)*value(i-1,:)' ; 
        % problem when the action is 'do nothing'
        minp_3=cxu(j,'maintenance')+alpha*transitionMatrices(3).matrix(j+1,:)*value(i-1,:)';
        end
        
        % Formulating Bellman equation and assigning the solution to value
        % matrix
        [value(i,j+1),idx]=min([minp_1,minp_2,minp_3]); %V(i,j)=min(...)
        
        if idx==1
            decision(i,j+1)={'replacement'};
        elseif idx==2
            decision(i,j+1)={'do nothing'};
        elseif idx==3
            decision(i,j+1)={'maintenance'};
        end
        
        if i>=2
            %Stores the difference between each entry of v(i+1) and v(i),
            %it acts like a history of differences.
            diffmat(i+2,:)=value(i,:)-value(i-1,:); 
        elseif i==1
            diffmat(i+1,:)=inf;
        end    
    end
    i=i+1;
end

solution.value=value(end,:);
solution.decision=decision(end,:);
%You code ends here
end