function [ time,POP] =Two_Arch2( problem_name )
% Usage: [ time,POP] =Two_Arch2( problem_name )
%
% Input:
% problem_name  - Benchmark MOP
%
% Output: 
% time          - Execution Time
% POP           - Obtained Solution Set with c Decision Variables and m Objectives
%
    %%%%    Authors:    Handing Wang, Licheng Jiao, Xin Yao
    %%%%    Xidian University, China, and University of Birmingham, UK
    %%%%    EMAIL:      wanghanding.patch@gmail.com, X.Yao@cs.bham.ac.uk
    %%%%    WEBSITE:    http://www.cs.bham.ac.uk/~xin/
    %%%%    DATE:       August 2014
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Licheng Jiao, Xin Yao, An Improved Two-Archive Algorithm for Many-Objective Optimization, Evolutionary Computation, IEEE Transactions on, Accepted, 10.1109/TEVC.2014.2350987.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------
rand('state',sum(100*clock));
%Problem Settings
%gmax- Maximum Iterations
%c-No. of Decision Variables
%n-Scale for CA and DA
switch problem_name
   case {'ZDT1','ZDT2','ZDT3'}
        gmax=100;
        c=30;
        n=100;
    case {'ZDT6','ZDT4'}
        gmax=100;
        c=10;
        n=100;
    case {'DTLZ2_3','DTLZ2','DTLZ4'}
        gmax=100;
        c=12;
        n=100;
    case {'DTLZ1','DTLZ3'}
        gmax=200;
        c=7;
        n=100;
    case {'DTLZ1_2','DTLZ3_2'}
        gmax=300;
        c=6;
        n=100;
    case {'DTLZ1_4','DTLZ3_4'}
        gmax=300;
        c=8;
        n=100;
    case {'DTLZ1_5','DTLZ3_5'}
        gmax=300;
        c=9;
        n=100;
    case {'DTLZ1_6','DTLZ3_6'}
        gmax=300;
        c=10;
        n=100;
    case {'DTLZ1_7','DTLZ3_7'}
        gmax=300;
        c=11;
        n=100;
    case {'DTLZ1_8','DTLZ3_8'}
        gmax=300;
        c=12;
        n=100;
    case {'DTLZ1_9','DTLZ3_9'}
        gmax=300;
        c=13;
        n=100;
    case {'DTLZ1_10','DTLZ3_10'}
        gmax=300;
        c=14;
        n=100;
    case {'DTLZ1_15','DTLZ3_15'}
        gmax=300;
        c=19;
        n=200;
    case {'DTLZ1_20','DTLZ3_20'}
        gmax=300;
        c=24;
        n=200;    
    case {'DTLZ2_2','DTLZ4_2'}
        gmax=300;
        c=11;
        n=100;
    case {'DTLZ2_4','DTLZ4_4'}
        gmax=300;
        c=13;
        n=100;
    case {'DTLZ2_5','DTLZ4_5'}
        gmax=300;
        c=14;
        n=100;
    case {'DTLZ2_6','DTLZ4_6'}
        gmax=300;
        c=15;
        n=100;
    case {'DTLZ2_7','DTLZ4_7'}
        gmax=300;
        c=16;
        n=100;
    case {'DTLZ2_8','DTLZ4_8'}
        gmax=300;
        c=17;
        n=100;
    case {'DTLZ2_9','DTLZ4_9'}
        gmax=300;
        c=18;
        n=100;
    case {'DTLZ2_10','DTLZ4_10'}
        gmax=300;
        c=19;
        n=100;
    case {'DTLZ2_15','DTLZ4_15'}
        gmax=300;
        c=24;
        n=200;
    case {'DTLZ2_20','DTLZ4_20'}
        gmax=300;
        c=29;
        n=200;
    case {'WFG1_2','WFG1_3','WFG1_4','WFG1_5','WFG1_6','WFG1_7','WFG1_8','WFG1_9','WFG1_10'}
        gmax=300;
        c=20;
        n=100;
    case {'WFG2_2','WFG2_3','WFG2_4','WFG2_5','WFG2_6','WFG2_7','WFG2_8','WFG2_9','WFG2_10'}
        gmax=300;
        c=20;
        n=100;
    case {'WFG3_2','WFG3_3','WFG3_4','WFG3_5','WFG3_6','WFG3_7','WFG3_8','WFG3_9','WFG3_10'}
        gmax=300;
        c=20;
        n=100;
    case {'WFG4_2','WFG4_3','WFG4_4','WFG4_5','WFG4_6','WFG4_7','WFG4_8','WFG4_9','WFG4_10'}
        gmax=300;
        c=20;
        n=100;
    case {'WFG5_2','WFG5_3','WFG5_4','WFG5_5','WFG5_6','WFG5_7','WFG5_8','WFG5_9','WFG5_10'}
        gmax=300;
        c=20;
        n=100;
    case {'WFG6_2','WFG6_3','WFG6_4','WFG6_5','WFG6_6','WFG6_7','WFG6_8','WFG6_9','WFG6_10'}
        gmax=300;
        c=20;
        n=100;
    case {'WFG7_2','WFG7_3','WFG7_4','WFG7_5','WFG7_6','WFG7_7','WFG7_8','WFG7_9','WFG7_10'}
        gmax=300;
        c=20;
        n=100;
    case {'WFG8_2','WFG8_3','WFG8_4','WFG8_5','WFG8_6','WFG8_7','WFG8_8','WFG8_9','WFG8_10'}
        gmax=300;
        c=20;
        n=100;
    case {'WFG9_2','WFG9_3','WFG9_4','WFG9_5','WFG9_6','WFG9_7','WFG9_8','WFG9_9','WFG9_10'}
        gmax=300;
        c=20;
        n=100;
end

pc=1;%Crossover Probability
pm=0.1;%Mutation Probability
[ bu,bd ] = generate_boundary( problem_name,c );%Upper and Lower Bound of Decision Variables
tic;
%Initialization
POP = initialize_pop(n,c,bu,bd);
obj =compute_objectives(POP,c,problem_name);
m=size(obj,2);%No. of Objectives
p=1/m;%Lp-norm
POP=[POP,obj];
g=1;%Iteration Counter
%Initial DA and CA
DA = find_nondominated( POP,c,m);
CA=[];
while g<=gmax
    %Crossover
    if size(CA,1)==0
        NPOP1=SBX( POP(:,1:c+m),bu,bd,pc,n );
        obj =compute_objectives(NPOP1,c,problem_name);
        NPOP1=[NPOP1,obj];
    else
        NPOP1 = SBXCD( CA,DA,bu,bd,pc,n,m  );
        obj =compute_objectives(NPOP1,c,problem_name);
        NPOP1=[NPOP1,obj];
    end
    %Mutation
    if size(CA,1)==0
        NPOP2=mutation(DA(:,1:c+m),bu,bd,pm,n);
        obj =compute_objectives(NPOP2,c,problem_name);
        NPOP2=[NPOP2,obj];
    else
        NPOP2=mutation(CA,bu,bd,pm,n);
        obj =compute_objectives(NPOP2,c,problem_name);
        NPOP2=[NPOP2,obj];
    end

    NPOP=[NPOP1;NPOP2];
    CPOP=find_nondominated( NPOP,c,m);
    %Update Both CA and DA
    [ CA,DA] = updateTA(CPOP,CA,NPOP,DA,c,m,n,p);
    g=g+1;
end

POP=DA;%Final Output
toc;
time=toc;
end

