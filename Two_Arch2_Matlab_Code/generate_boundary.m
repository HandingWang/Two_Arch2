function [ bu,bd ] = generate_boundary( problem_name,c )
% Usage: [ bu,bd ] = generate_boundary( problem_name,c )
%
% Input:
% problem_name  - Benchmark MOP (ZDT, DTLZ, WFG Problems)
% c             -No. of Decision Variables
%
% Output: 
% bu            - Upper Bound
% bd            - Lower Bound
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
switch problem_name
    case {'DTLZ1','DTLZ2','DTLZ3','DTLZ4','DTLZ1_2','DTLZ1_4','DTLZ1_5','DTLZ1_6','DTLZ1_7','DTLZ1_8','DTLZ1_9','DTLZ1_10','DTLZ1_15','DTLZ2_2','DTLZ2_3','DTLZ2_4','DTLZ2_5','DTLZ2_6','DTLZ2_7','DTLZ2_8','DTLZ2_9','DTLZ2_10','DTLZ2_15','DTLZ2_20','DTLZ3_2','DTLZ3_3','DTLZ3_4','DTLZ3_5','DTLZ3_6','DTLZ3_7','DTLZ3_8','DTLZ3_9','DTLZ3_10','DTLZ4_2','DTLZ4_3','DTLZ4_4','DTLZ4_5','DTLZ4_6','DTLZ4_7','DTLZ4_8','DTLZ4_9','DTLZ4_10','DTLZ3_15','DTLZ4_15','DTLZ4_20','DTLZN_2','DTLZN_3','DTLZN_4','DTLZN_5','DTLZN_6','DTLZN_7','DTLZN_8','DTLZN_9','DTLZN_10','DTLZN_15','DTLZN_20','DTLZ1_20','DTLZ3_20'}
        bu=ones(1,c);
        bd=zeros(1,c);
    case 'ZDT4',
        bu(1)=1;
        bd(1)=0;
        bu(2:c)=5*ones(1,c-1);
        bd(2:c)=-5*ones(1,c-1);
    case {'WFG1_2','WFG1_3','WFG1_4','WFG1_5','WFG1_6','WFG1_7','WFG1_8','WFG1_9','WFG1_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
    case {'WFG2_2','WFG2_3','WFG2_4','WFG2_5','WFG2_6','WFG2_7','WFG2_8','WFG2_9','WFG2_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
    case {'WFG3_2','WFG3_3','WFG3_4','WFG3_5','WFG3_6','WFG3_7','WFG3_8','WFG3_9','WFG3_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
    case {'WFG4_2','WFG4_3','WFG4_4','WFG4_5','WFG4_6','WFG4_7','WFG4_8','WFG4_9','WFG4_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
    case {'WFG5_2','WFG5_3','WFG5_4','WFG5_5','WFG5_6','WFG5_7','WFG5_8','WFG5_9','WFG5_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
    case {'WFG6_2','WFG6_3','WFG6_4','WFG6_5','WFG6_6','WFG6_7','WFG6_8','WFG6_9','WFG6_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
    case {'WFG7_2','WFG7_3','WFG7_4','WFG7_5','WFG7_6','WFG7_7','WFG7_8','WFG7_9','WFG7_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
    case {'WFG8_2','WFG8_3','WFG8_4','WFG8_5','WFG8_6','WFG8_7','WFG8_8','WFG8_9','WFG8_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
    case {'WFG9_2','WFG9_3','WFG9_4','WFG9_5','WFG9_6','WFG9_7','WFG9_8','WFG9_9','WFG9_10'}
        bd=zeros(1,c);       
        bu=2:2:2*c;
end
end

