function [ POP ] = initialize_pop(n,c,bu,bd)
% Usage: [ POP ] = initialize_pop(n,c,bu,bd)
%
% Input:
% bu            -Upper Bound
% bd            -Lower Bound
% c             -No. of Decision Variables
% n             -Population Scale
%
% Output: 
% POP           -Initial Population
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
POP=rand(n,c).*(ones(n,1)*(bu-bd))+ones(n,1)*bd;
end