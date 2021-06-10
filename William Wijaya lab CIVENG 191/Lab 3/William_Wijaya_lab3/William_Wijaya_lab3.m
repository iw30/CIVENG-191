clear all;clc;close all
for i=1:6
    if i==1
        load('graph1.mat')
        visGraph(graph.edges)
    elseif i==2
        load('graph2.mat')
        visGraph(graph.edges)
    elseif i==3
        load('graph3.mat')
        visGraph(graph.edges)
    elseif i==4
        load('graph4.mat')
        visGraph(graph.edges)
    elseif i==5
        load('graph5.mat')
        visGraph(graph.edges)
    else 
        load('graph6.mat')
        visGraph(graph.edges)
    end


prob=formLP(graph.edges)

prob=solveLP(prob)

if prob.isFeasible==1
    visGraph([prob.sol_edges,prob.cost_edges])
end
    
fprintf('Press enter to continue next problem:\n');pause;
clc;clear all;close all
end
fprintf('Checkpoint 1 lab 3 Done!\n')