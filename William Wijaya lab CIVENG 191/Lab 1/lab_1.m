
% William Wijaya/CIVENG 191/lab 1
% Run this matlab file to start the program

%% 2

% Data

clear all,clc;

% Maximum hardness data
h1=250;  %Source 1
h2=200;  %Source 2
h3=2300; %Source 3
h4=700;  %Source 3

% Daily water need
n1=30;   %Appletown
n2=10;   %Berrytown
n3=50;   %Cherrytown
n4=20;   %Grapetown
n5=40;   %Mangotown

% Supply limit per day
s1=15;   %Source1
s2=10;   %Source2
s3=55;   %Source3
s4=85;   %Source4

% Matrix
O=ones(1,5);  %ones 1x5 matrix
Z=zeros(1,5); %zeros 1x5 matrix

% lower bound
lb=zeros(20,1);

% Constructing "c^T" matrix

%Cost of providing water from source 1
c1=[365,375,355,360,370];
%Cost of providing water from source 2
c2=[425,430,435,440,420];
%Cost of providing water from source 3
c3=[705,720,690,725,730];
%Cost of providing water from source 4
c4=[2005,2020,2000,1990,2025];
%Concatenate matrix to c^T
c=[c1,c2,c3,c4];

% Constructing "A" Matrix

% Daily need of Appletown
T_A=[1,0,0,0,0];
N1=[T_A,T_A,T_A,T_A];
% Daily need of Berrytown
T_B=[0,1,0,0,0];
N2=[T_B,T_B,T_B,T_B];
% Daily need of Cherrytown
T_C=[0,0,1,0,0];
N3=[T_C,T_C,T_C,T_C];
% Daily need of Grapetown
T_G=[0,0,0,1,0];
N4=[T_G,T_G,T_G,T_G];
% Daily need of Mangotown
T_M=[0,0,0,0,1];
N5=[T_M,T_M,T_M,T_M];
% Concatenate Daily need constraint
N=[N1;N2;N3;N4;N5;-1.*N1;-1.*N2;-1.*N3;-1.*N4;-1.*N5];

% Supply limit of Source 1
S1=[O,Z,Z,Z];
% Supply limit of Source 2
S2=[Z,O,Z,Z];
% Supply limit of Source 3
S3=[Z,Z,O,Z];
% Supply limit of Source 4
S4=[Z,Z,Z,O];
% Concatenate supply constraint
S=[S1;S2;S3;S4];

% Common pipe constraint
P=[O,O,Z,Z];

% Maximum Hardness for Appletown
H1=[h1.*T_A,h2.*T_A,h3.*T_A,h4.*T_A];
% Maximum Hardness for Berrytown
H2=[h1.*T_B,h2.*T_B,h3.*T_B,h4.*T_B];
% Maximum Hardness for Cherrytown
H3=[h1.*T_C,h2.*T_C,h3.*T_C,h4.*T_C];
% Maximum Hardness for Grapetown
H4=[h1.*T_G,h2.*T_G,h3.*T_G,h4.*T_G];
% Maximum Hardness for Mangotown
H5=[h1.*T_M,h2.*T_M,h3.*T_M,h4.*T_M];
% Concatenate Maximum Hardness constraint
H=[H1;H2;H3;H4;H5];

% Concatenate "A" Matrix
A=[N;S;P;H];

% Constructing "b" Matrix

% Daily need matrix
n=[n1,n2,n3,n4,n5];

% Supply limit matrix
s=[s1,s2,s3,s4];

% Maximum Hardness matrix
h=1200.*n;

% Concatenate to "b" matrix
b=[n,-1.*n,s,20,h];

% Solve Linear Programming
options = optimoptions('linprog','Algorithm','dual-simplex');
[x,fval,exitflag,output,lambda]=linprog(c,A,b,[],[],lb,[],options);

% Constructing a table summarizing the answer

nameS1=["Source 1 to Appletown","Source 1 to Berrytown",...
    "Source 1 to Cherrytown","Source 1 to Grapetown",...
    "Source 1 to Mangotown"];
nameS2=["Source 2 to Appletown","Source 2 to Berrytown",...
    "Source 2 to Cherrytown","Source 2 to Grapetown",...
    "Source 2 to Mangotown"];
nameS3=["Source 3 to Appletown","Source 3 to Berrytown",...
    "Source 3 to Cherrytown","Source 3 to Grapetown",...
    "Source 3 to Mangotown"];
nameS4=["Source 4 to Appletown","Source 4 to Berrytown",...
    "Source 4 to Cherrytown","Source 4 to Grapetown",...
    "Source 4 to Mangotown"];
Deliver_route=[nameS1,nameS2,nameS3,nameS4]';
Amount_of_water_delivered=x;
fprintf("The result for no 2 is displayed below in a table:\n")
result_2=table(Deliver_route,Amount_of_water_delivered)
mincost_2=fval;
fprintf("The minimum transportation cost is $%4.2f \n\n",mincost_2)

%Checking for active constraints
fprintf("The table for constraints is shown below: \n") 
activeValue=A*x-b'; 
nameC1=["Daily need of Appletown (+)","Daily need of Berrytown (+)",...
    "Daily need of Cherrytown (+)","Daily need of Grapetown (+)",...
    "Daily need of Mangotown (+)"];
nameC2=["Daily need of Appletown (-)","Daily need of Berrytown (-)",...
    "Daily need of Cherrytown (-)","Daily need of Grapetown (-)",...
    "Daily need of Mangotown (-)"];
nameC3=["Supply limit of S1","Supply limit of S2","Supply limit of S3",...
    "Supply limit of S4"];
nameC4=["Common pipe"];
nameC5=["Maximum Hardness of Appletown","Maximum Hardness of Berrytown",...
    "Maximum Hardness of Cherrytown","Maximum Hardness of Grapetown",...
    "Maximum Hardness of Mangotown"];
constraint=[nameC1,nameC2,nameC3,nameC4,nameC5]';
active_or_inactive=["active","active","active","active","active",...
    "active","active","active","active","active",...
    "active","inactive","inactive","inactive","active",...
    "active","active","active","active","active"]';
tableConstraint=table(constraint,activeValue,active_or_inactive)

%% 3a

display("Press enter to continue no 3a:");pause;clc;

%Decrease the demand of water by 5% from each town
n=0.95.*n;

b=[n,-1.*n,s,20,h];
options = optimoptions('linprog','Algorithm','dual-simplex');
[x,fval,exitflag,output,lambda]=linprog(c,A,b,[],[],lb,[],options);
Amount_of_water_delivered=x;
fprintf("The result for no 2 is displayed below in a table:\n")
result_3a=table(Deliver_route,Amount_of_water_delivered)
fprintf("The minimum transportation cost is $%4.2f \n",fval)

%% 3b
display("Press enter to continue no 3b:");pause;clc;

% Increase the demand of water by 5% from each town
n=[n1,n2,n3,n4,n5];
n=1.05.*n;

b=[n,-1.*n,s,20,h];
options = optimoptions('linprog','Algorithm','dual-simplex');
[x,fval,exitflag,output,lambda]=linprog(c,A,b,[],[],lb,[],options)
fprintf("The result for no 3b is not feasible,")
fprintf(" i.e., no optimal sol found \n")

%% 3c
display("Press enter to continue no 3c:");pause;clc;
n=[n1,n2,n3,n4,n5];
s=[s1,s2,95,75];
b=[n,-1.*n,s,20,h];
options = optimoptions('linprog','Algorithm','dual-simplex');
[x,fval,exitflag,output,lambda]=linprog(c,A,b,[],[],lb,[],options)
fprintf("The result for no 3c is not feasible,")
fprintf(" i.e., no optimal sol found \n")

%% 4

% First, we note that sending water to Kiwi is the same as sending it to 
% Berrytown. So, we just need to know what is the price/ML to send water to
% Berrytown. Now, from the table number 2, we can extract the amount of
% water delivered from source i to Berrytown (i=1,2,3,4)

display("Press enter to continue no 4:");pause;clc;
s=[s1,s2,s3,s4];
b=[n,-1.*n,s,20,h];
options = optimoptions('linprog','Algorithm','dual-simplex');
[x,fval,exitflag,output,lambda]=linprog(c,A,b,[],[],lb,[],options);

% Cost of sending water to Berrytown per ML
cost_B=(x(2)*c(2)+x(7)*c(7)+x(12)*c(12)+x(17)*c(17))/n2;
fprintf("Cost of sending water to Berrytown per ML is $%4.2f \n",cost_B)

if cost_B<1500
    fprintf("Accept Kiwi Inc's offer since this is a profit ")
    fprintf("(price/ML<1500) \n")
else 
    fprintf("Reject Kiwi Inc's offer since this is a loss ")
    fprintf("(price/ML>1500) \n")
end

%% 5

display("Press enter to continue no 5:");pause;clc;
% The program is the same as no 2, but the only changes are:
% I increase s1 by 2
% And, I added the mincost_5 by (5000/365)

s1=s1+2;
s=[s1,s2,s3,s4];
b=[n,-1.*n,s,20,h];

% Solve Linear Programming
options = optimoptions('linprog','Algorithm','dual-simplex');
[x,fval,exitflag,output,lambda]=linprog(c,A,b,[],[],lb,[],options);

% Constructing a table summarizing the answer

Amount_of_water_delivered=x;
fprintf("The result for no 5 is displayed below in a table: \n")
result_5=table(Deliver_route,Amount_of_water_delivered)
mincost_5=fval+(5000/365);
fprintf("The minimum transportation cost is $%4.2f \n",mincost_5)
difference=mincost_2 - mincost_5;
fprintf("By adding a small structure to source 1,")
fprintf(" we gain profit of $%4.2f \n",difference)
fprintf("Choose this option to gain the profit \n")

%% 6
display("Press enter to continue no 6:");pause;clc;

% Maximum hardness data
h1=250;  %Source 1
h2=200;  %Source 2
h3=2300; %Source 3
h4=700;  %Source 3

% Daily water need
n1=30;   %Appletown
n2=10;   %Berrytown
n3=50;   %Cherrytown
n4=20;   %Grapetown
n5=40;   %Mangotown

% Supply limit per day
s1=15;   %Source1
s2=10;   %Source2
s3=56;   %Source3, add this by one
s4=85;   %Source4

s=[s1,s2,s3,s4];
n=[n1,n2,n3,n4,n5];
b=[n,-1.*n,s,20,h];

options = optimoptions('linprog','Algorithm','dual-simplex');
[x,fval,exitflag,output,lambda]=linprog(c,A,b,[],[],lb,[],options);

Amount_of_water_delivered=x;
fprintf("The result for no 6 is displayed below in a table:\n")
result_6=table(Deliver_route,Amount_of_water_delivered)
mincost_6=fval;
fprintf("The minimum transportation cost for no 6 is $%4.2f \n",mincost_6)
fprintf("The minimum transportation cost for no 2 is $%4.2f \n",mincost_2)
fprintf("The minimum cost for no 2 and 6 are the same,")
different= mincost_6-mincost_2;
fprintf(" with a difference of %4.2f \n", different)
fprintf("So, do not pay for any increase of supply limit 3 because the")
fprintf(" minimum cost will not change \n\n")

fprintf("Lab 1 finished! \n")