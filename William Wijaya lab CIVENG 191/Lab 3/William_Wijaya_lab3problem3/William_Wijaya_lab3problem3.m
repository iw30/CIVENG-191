clear all;clc;close all
%This is for problem 3 only!
for i=1:6
    if i==1
        load('graph1.mat')
    elseif i==2
        load('graph2.mat')
    elseif i==3
        load('graph3.mat')
    elseif i==4
        load('graph4.mat')
    elseif i==5
        load('graph5.mat')
    else 
        load('graph6.mat')
    end

[tour,cost]=solveTSP(graph.edges);
fprintf('\nPress enter to continue next problem:\n');pause;
clc;clear all;close all
end
fprintf('lab 3 problem 3 Done!\n')