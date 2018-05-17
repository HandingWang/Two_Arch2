function [ NPOP ] = find_nondominated( POP,c,m)
% Usage: [ NPOP ] = find_nondominated( POP,c,m)
%
% Input:
% m             -No. of Objectives
% c             -No. of Decision Variables
% POP           -Population of Decision Variables
%
% Output: 
% NPOP          - Non-dominated solution set
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
n=size(POP,1);
i=1;
while i<=size(POP,1)
    flag=0;
    j=i+1;
    while j<=size(POP,1)
        x=dominated_relationship(POP(i,:),POP(j,:),m,c);
        if x==2
            flag=1;
            break;
        elseif x==3
            POP(j,:)=[];
        elseif x==1
            POP(j,:)=[];
        else
            j=j+1;
        end
    end
    if flag==1
        POP(i,:)=[];
    else
        i=i+1;
    end
end
NPOP=POP;
end



