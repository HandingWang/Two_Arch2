function [x] = dominated_relationship(a,b,m,c)
% Usage: [x] = dominated_relationship(a,b,m,c)
%
% Input:
% a             - Solution a
% b             - Solution b
% m             - No. of Objectives
% c             - No. of Decision Variables
%
% Output: 
% x             - dominance relationship (1-a dominates b, 2-b dominates a, 3-a = b, 4-a and b are non-dominated by each other)
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
t=0;
q=0;
p=0;
for i=1:m
    if a(1,c+i)<=b(1,c+i)
        t=t+1;
    end
    if  a(1,c+i)>= b(1,c+i)
        q=q+1;
    end
    if  a(1,c+i)== b(1,c+i)
        p=p+1;
    end
end

if t==m&p~=m
    x=1;
elseif q==m&p~=m
    x=2;
elseif p==m
    x=3;
else
    x=4;
end
end

