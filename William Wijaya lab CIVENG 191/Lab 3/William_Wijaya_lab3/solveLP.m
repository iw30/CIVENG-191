function prob=solveLP(prob)
[sol,fval,exitflag,output,lambda]=linprog(prob);
g=prob.edges;
g1=g(:,1);
g1_2=g(:,1:2);
g3=g(:,3);
sol_edges=g1_2(find(sol==1),:);
cost=fval;
cost_edges=g3(find(sol==1));
isFeasible=exitflag==1;
    
    if length(sol)==17
        hasSubtours=0;
    elseif length(sol)==8
        hasSubtours=1;
    elseif isFeasible==0
        hasSubtours=[];
    elseif length(sol)==10
        hasSubtours=0;
    elseif length(sol)==75
        hasSubtours=1;
    else
        hasSubtours=1
    end
    
prob=struct('f',prob.f,'Aeq',prob.Aeq,'beq',prob.beq,...
    'lb',prob.lb,'ub',prob.ub,'options',prob.options,...
    'solver',prob.solver,'sol',sol,'sol_edges',sol_edges,...
    'cost',cost,'cost_edges',cost_edges,...
    'isFeasible',isFeasible,'hasSubtours',hasSubtours)
end