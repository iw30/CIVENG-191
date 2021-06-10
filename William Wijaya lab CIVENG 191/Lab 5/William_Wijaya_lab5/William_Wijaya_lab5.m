%% This section is for problem 1 in lab 5
%% Insert pq command
clear all,clc,close all

a=0:3;
n=10.^a;


for R=1:500 %Repeat the process "R" times
    for i=0:3 % for n=10^i
    pq=PriorityQueue();
    % Inserting "n" element in the priority queue starting from priority 2
        for j=2:n(i+1)+1 
        pq.Insert(j-1,j);
        end
        start=tic; % Start timer
        pq.Insert(j,1); %Inserting a first priority element into priority queue
        stop=toc(start); % Stop timer 
        tEnd_Insert(i+1,R)=stop; % Record time
        %clear pq % Clear priority queue and start from the beginning
     end
end

tEnd_Insert=mean(tEnd_Insert,2);
semilogx(n,tEnd_Insert,'r-x')
xlabel('Number of elements, n')
ylabel('Time')
title('Insert()')

%% ExtractMin command
clear all,clc,close all
a=0:3;
n=10.^a;

for R=1:500
    for i=0:3
        pq=PriorityQueue();
        for j=1:n(i+1)
            pq.Insert(j,j);
        end
        start=tic;
        pq.ExtractMin();
        stop=toc(start);
        tEnd_Insert(i+1,R)=stop;
    end
end

tEnd_Insert=mean(tEnd_Insert,2)
semilogx(n,tEnd_Insert,'r-x')
xlabel('Number of elements, n')
ylabel('Time')
title('ExtractMin()')

%% Decrease key command
clear all,clc,close all
a=0:3;
n=10.^a;

for R=1:500
    for i=0:3
        pq=PriorityQueue;
            for j=2:n(i+1)+1
            pq.Insert(j-1,j);
            end
        start=tic;
        pq.DecreaseKey(j-1,1);
        stop=toc(start);
        tEnd_Insert(i+1,R)=stop;
        clear pq
    end
end

tEnd_Insert=mean(tEnd_Insert,2)
semilogx(n,tEnd_Insert,'r-x')
xlabel('Number of elements, n')
ylabel('Time')
title('DecreaseKey()')

fprintf("\nDone!\n")