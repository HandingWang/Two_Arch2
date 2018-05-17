function [ obj ] =compute_objectives(POP,c,problem_name)
% Usage: [ obj ] =compute_objectives(POP,c,problem_name)
%
% Input:
% problem_name  - Benchmark MOP (ZDT, DTLZ, WFG Problems)
% c             -No. of Decision Variables
% POP           -Population of Decision Variables
%
% Output: 
% obj           - Calculated Objective Values
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
obj=[];
n=size(POP,1);
switch problem_name
    case 'ZDT1',
        obj(:,1)=POP(:,1);
        g=ones(n,1)+9*sum(POP(:,2:c),2)/(c-1);
        obj(:,2)=g.*(ones(n,1)-(POP(:,1)./g).^0.5);
    case 'ZDT2',
        obj(:,1)=POP(:,1);
        g=ones(n,1)+9*sum(POP(:,2:c),2)/(c-1);
        obj(:,2)=g.*(ones(n,1)-(POP(:,1)./g).^2);
    case 'ZDT3',
        obj(:,1)=POP(:,1);
        g=ones(n,1)+9*sum(POP(:,2:c),2)/(c-1);
        obj(:,2)=g.*(ones(n,1)-(POP(:,1)./g).^0.5-(POP(:,1)./g).*(sin(10*pi*POP(:,1))));
    case 'ZDT4',
        obj(:,1)=POP(:,1);
        g=ones(n,1)+10*(c-1)*ones(n,1)+sum((POP(:,2:c).^2-10*cos(4*pi*POP(:,2:c))),2);
        obj(:,2)=g.*(ones(n,1)-(POP(:,1)./g).^0.5);
    case 'ZDT6',
        obj(:,1)=ones(n,1)-exp(-4*POP(:,1)).*((sin(6*pi*POP(:,1))).^6);
        g=ones(n,1)+9*((sum(POP(:,2:c),2)/(c-1)).^0.25);
        obj(:,2)=g.*(ones(n,1)-obj(:,1).^2);
    case 'DTLZ1'
        g=100*((c-2)*ones(n,1)+sum((POP(:,3:c)-0.5).^2-cos(20*pi*(POP(:,3:c)-0.5)),2));
        obj(:,1)=POP(:,1).*POP(:,2).*(1+g);
        obj(:,2)=POP(:,1).*(1-POP(:,2)).*(1+g);
        obj(:,3)=(1-POP(:,1)).*(1+g);
    case 'DTLZ2'
        g=sum((POP(:,3:c)-0.5).^2,2);
        obj(:,1)=cos(0.5*pi*POP(:,1)).*cos(0.5*pi*POP(:,2)).*(1+g);
        obj(:,2)=cos(0.5*pi*POP(:,1)).*sin(0.5*pi*POP(:,2)).*(1+g);
        obj(:,3)=sin(0.5*pi*POP(:,1)).*(1+g);
    case 'DTLZ3'
        g=100*((c-2)*ones(n,1)+sum((POP(:,3:c)-0.5).^2-cos(20*pi*(POP(:,3:c)-0.5)),2));
        obj(:,1)=cos(0.5*pi*POP(:,1)).*cos(0.5*pi*POP(:,2)).*(1+g);
        obj(:,2)=cos(0.5*pi*POP(:,1)).*sin(0.5*pi*POP(:,2)).*(1+g);
        obj(:,3)=sin(0.5*pi*POP(:,1)).*(1+g);
    case 'DTLZ4'
        g=sum((POP(:,3:c)-0.5).^2,2);
        obj(:,1)=cos(0.5*pi*POP(:,1).^100).*cos(0.5*pi*POP(:,2).^100).*(1+g);
        obj(:,2)=cos(0.5*pi*POP(:,1).^100).*sin(0.5*pi*POP(:,2).^100).*(1+g);
        obj(:,3)=sin(0.5*pi*POP(:,1).^100).*(1+g);
    case 'DTLZ1_2'
        m=2;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_4'
        m=4;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_5'
        m=5;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_6'
        m=6;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_7'
        m=7;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_8'
        m=8;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_9'
        m=9;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_10'
        m=10;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_15'
        m=15;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);
    case 'DTLZ1_20'
        m=20;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(POP(:,1:m-1),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(POP(:,1:m-i),2).*(1-POP(:,m+1-i)).*(1+g);
        end
        obj(:,m)=(1-POP(:,1)).*(1+g);    
    case 'DTLZ2_2'
        m=2;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ2_3'
        m=3;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ2_4'
        m=4;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ2_5'
        m=5;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ2_6'
        m=6;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ2_7'
        m=7;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ2_8'
        m=8;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ2_9'
        m=9;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ2_10'
        m=10;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ2_15'
        m=15;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ2_20'
        m=20;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ3_2'
        m=2;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_3'
        m=3;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_4'
        m=4;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_5'
        m=5;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_6'
        m=6;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_7'
        m=7;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_8'
        m=8;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_9'
        m=9;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_10'
        m=10;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
     case 'DTLZ3_15'
        m=15;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);
    case 'DTLZ3_20'
        m=20;
        g=100*((c-m+1)*ones(n,1)+sum((POP(:,m:c)-0.5).^2-cos(20*pi*(POP(:,m:c)-0.5)),2));
        obj(:,1)=prod(cos(0.5*pi*POP(:,1:m-1)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod(cos(0.5*pi*POP(:,1:m-i)),2).*(sin(0.5*pi*POP(:,m+1-i))).*(1+g);
        end
        obj(:,m)=(sin(0.5*pi*POP(:,1))).*(1+g);    
     case 'DTLZ4_2'
        m=2;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_3'
        m=3;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_4'
        m=4;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_5'
        m=5;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_6'
        m=6;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_7'
        m=7;
       g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_8'
        m=8;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_9'
        m=9;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_10'
        m=10;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_15'
        m=15;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
     case 'DTLZ4_20'
        m=20;
        g=sum((POP(:,m:c)-0.5).^2,2);
        obj(:,1)=prod((cos(0.5*pi*POP(:,1:m-1).^100)),2).*(1+g);
        for i=2:m-1
            obj(:,i)=prod((cos(0.5*pi*POP(:,1:m-i).^100)),2).*((sin(0.5*pi*POP(:,m+1-i).^100))).*(1+g);
        end
        obj(:,m)=((sin(0.5*pi*POP(:,1).^100))).*(1+g);
    case 'WFG1_2'
        obj = WFG1(2,10,POP );
    case 'WFG1_3'
        obj = WFG1(3,10,POP );
    case 'WFG1_4'
        obj = WFG1(4,9,POP );
    case 'WFG1_5'
        obj = WFG1(5,12,POP );
    case 'WFG1_6'
        obj = WFG1(6,10,POP );
    case 'WFG1_7'
        obj = WFG1(7,12,POP );
    case 'WFG1_8'
        obj = WFG1(8,7,POP );
    case 'WFG1_9'
        obj = WFG1(9,8,POP );
    case 'WFG1_10'
        obj = WFG1(10,9,POP );
    case 'WFG2_2'
        obj = WFG2(2,10,POP );
    case 'WFG2_3'
        obj = WFG2(3,10,POP );
    case 'WFG2_4'
        obj = WFG2(4,9,POP );
    case 'WFG2_5'
        obj = WFG2(5,12,POP );
    case 'WFG2_6'
        obj = WFG2(6,10,POP );
    case 'WFG2_7'
        obj = WFG2(7,12,POP );
    case 'WFG2_8'
        obj = WFG2(8,7,POP );
    case 'WFG2_9'
        obj = WFG2(9,8,POP );
    case 'WFG2_10'
        obj = WFG2(10,9,POP );
    case 'WFG3_2'
        obj = WFG3(2,10,POP );
    case 'WFG3_3'
        obj = WFG3(3,10,POP );
    case 'WFG3_4'
        obj = WFG3(4,9,POP );
    case 'WFG3_5'
        obj = WFG3(5,12,POP );
    case 'WFG3_6'
        obj = WFG3(6,10,POP );
    case 'WFG3_7'
        obj = WFG3(7,12,POP );
    case 'WFG3_8'
        obj = WFG3(8,7,POP );
    case 'WFG3_9'
        obj = WFG3(9,8,POP );
    case 'WFG3_10'
        obj = WFG3(10,9,POP );
    case 'WFG4_2'
        obj = WFG4(2,10,POP );
    case 'WFG4_3'
        obj = WFG4(3,10,POP );
    case 'WFG4_4'
        obj = WFG4(4,9,POP );
    case 'WFG4_5'
        obj = WFG4(5,12,POP );
    case 'WFG4_6'
        obj = WFG4(6,10,POP );
    case 'WFG4_7'
        obj = WFG4(7,12,POP );
    case 'WFG4_8'
        obj = WFG4(8,7,POP );
    case 'WFG4_9'
        obj = WFG4(9,8,POP );
    case 'WFG4_10'
        obj = WFG4(10,9,POP );
    case 'WFG5_2'
        obj = WFG5(2,10,POP );
    case 'WFG5_3'
        obj = WFG5(3,10,POP );
    case 'WFG5_4'
        obj = WFG5(4,9,POP );
    case 'WFG5_5'
        obj = WFG5(5,12,POP );
    case 'WFG5_6'
        obj = WFG5(6,10,POP );
    case 'WFG5_7'
        obj = WFG5(7,12,POP );
    case 'WFG5_8'
        obj = WFG5(8,7,POP );
    case 'WFG5_9'
        obj = WFG5(9,8,POP );
    case 'WFG5_10'
        obj = WFG5(10,9,POP );
    case 'WFG6_2'
        obj = WFG6(2,10,POP );
    case 'WFG6_3'
        obj = WFG6(3,10,POP );
    case 'WFG6_4'
        obj = WFG6(4,9,POP );
    case 'WFG6_5'
        obj = WFG6(5,12,POP );
    case 'WFG6_6'
        obj = WFG6(6,10,POP );
    case 'WFG6_7'
        obj = WFG6(7,12,POP );
    case 'WFG6_8'
        obj = WFG6(8,7,POP );
    case 'WFG6_9'
        obj = WFG6(9,8,POP );
    case 'WFG6_10'
        obj = WFG6(10,9,POP );
    case 'WFG7_2'
        obj = WFG7(2,10,POP );
    case 'WFG7_3'
        obj = WFG7(3,10,POP );
    case 'WFG7_4'
        obj = WFG7(4,9,POP );
    case 'WFG7_5'
        obj = WFG7(5,12,POP );
    case 'WFG7_6'
        obj = WFG7(6,10,POP );
    case 'WFG7_7'
        obj = WFG7(7,12,POP );
    case 'WFG7_8'
        obj = WFG7(8,7,POP );
    case 'WFG7_9'
        obj = WFG7(9,8,POP );
    case 'WFG7_10'
        obj = WFG7(10,9,POP );
    case 'WFG8_2'
        obj = WFG8(2,10,POP );
    case 'WFG8_3'
        obj = WFG8(3,10,POP );
    case 'WFG8_4'
        obj = WFG8(4,9,POP );
    case 'WFG8_5'
        obj = WFG8(5,12,POP );
    case 'WFG8_6'
        obj = WFG8(6,10,POP );
    case 'WFG8_7'
        obj = WFG8(7,12,POP );
    case 'WFG8_8'
        obj = WFG8(8,7,POP );
    case 'WFG8_9'
        obj = WFG8(9,8,POP );
    case 'WFG8_10'
        obj = WFG8(10,9,POP );
    case 'WFG9_2'
        obj = WFG9(2,10,POP );
    case 'WFG9_3'
        obj = WFG9(3,10,POP );
    case 'WFG9_4'
        obj = WFG9(4,9,POP );
    case 'WFG9_5'
        obj = WFG9(5,12,POP );
    case 'WFG9_6'
        obj = WFG9(6,10,POP );
    case 'WFG9_7'
        obj = WFG9(7,12,POP );
    case 'WFG9_8'
        obj = WFG9(8,7,POP );
    case 'WFG9_9'
        obj = WFG9(9,8,POP );
    case 'WFG9_10'
        obj = WFG9(10,9,POP );
end


end

% WFG Problems
function [ obj ] =WFG1( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
t1(:,1:k)=x(:,1:k);
t1(:,k+1:c)=abs(x(:,k+1:c)-0.35)./abs(floor(0.35-x(:,k+1:c))+0.35);
t2(:,1:k)=t1(:,1:k);
A=0.8;
B=0.75;
C=0.85;
for i=k+1:c
    t2(:,i)=A+min(zeros(n,1),floor(t1(:,i)-B)).*(A*(B-t1(:,i))/B)-min(zeros(n,1),floor(C-t1(:,i))).*((1-A)*(t1(:,i)-C)/(1-C));
end
a=0.02;
for i=1:c
    t3(:,i)=t2(:,i).^a;
end
for i=1:m-1
    w=2*(i-1)*k/(m-1)+1:2:2*i*k/(m-1);
    y=t3(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    W=repmat(w,n,1);
    t4(:,i)=sum(y.*W,2)/sum(w,2);
end
w=2*(k+1):2:2*c;
y=t3(:,k+1:end);
W=repmat(w,n,1);
t4(:,m)=sum(y.*W,2)/sum(w,2);
x=[];
for i=1:1
    x(:,i)=max(t4(:,m),ones(n,1)).*(t4(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t4(:,m),ones(n,1)).*(t4(:,i)-0.5)+0.5;
end
x(:,m)=t4(:,m);

A=5;
a=1;

h(:,1)=prod(1-cos(0.5*pi*x(:,1:m-1)),2);
for i=2:m-1
    h(:,i)=prod(1-cos(0.5*pi*x(:,1:m-i)),2).*(1-sin(0.5*pi*x(:,m-i+1)));
end
h(:,m)=(1-x(:,1)-cos(2*A*pi*x(:,1)+0.5*pi)/(2*A*pi)).^a;
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end

function [ obj ] =WFG2( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
t1(:,1:k)=x(:,1:k);
t1(:,k+1:c)=abs(x(:,k+1:c)-0.35)./abs(floor(0.35-x(:,k+1:c))+0.35);
t2(:,1:k)=t1(:,1:k);
for i=k+1:k+(c-k)/2
    y=t1(:,k+(i-k)*2-1:k+(i-k)*2);
    t2(:,i)=(y(:,1)+y(:,2)+2*abs(y(:,1)-y(:,2)))/3;
end
for i=1:m-1
    y=t2(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    t3(:,i)=sum(y,2)/size(y,2);
end
y=t2(:,k+1:end);
t3(:,m)=sum(y,2)/size(y,2);
x=[];
for i=1:1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
x(:,m)=t3(:,m);
A=5;
a=1;
b=1;
h(:,1)=prod(1-cos(0.5*pi*x(:,1:m-1)),2);
for i=2:m-1
    h(:,i)=prod(1-cos(0.5*pi*x(:,1:m-i)),2).*(1-sin(0.5*pi*x(:,m-i+1)));
end
h(:,m)=1-(x(:,1).^a).*(cos(A*pi*(x(:,1).^b)).^2);
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end

function [ obj ] =WFG3( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
t1(:,1:k)=x(:,1:k);
t1(:,k+1:c)=abs(x(:,k+1:c)-0.35)./abs(floor(0.35-x(:,k+1:c))+0.35);
t2(:,1:k)=t1(:,1:k);
for i=k+1:k+(c-k)/2
    y=t1(:,k+(i-k)*2-1:k+(i-k)*2);
    t2(:,i)=(y(:,1)+y(:,2)+2*abs(y(:,1)-y(:,2)))/3;
end
for i=1:m-1
    y=t2(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    t3(:,i)=sum(y,2)/size(y,2);
end
y=t2(:,k+1:end);
t3(:,m)=sum(y,2)/size(y,2);
% t=max(t3(:,m),ones(n,1));
x=[];
for i=1:1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t3(:,m),zeros(n,1)).*(t3(:,i)-0.5)+0.5;
end
x(:,m)=t3(:,m);

h(:,1)=prod(x(:,1:m-1),2);
for i=2:m-1
    h(:,i)=prod(x(:,1:m-i),2).*(1-x(:,m-i+1));
end
h(:,m)=1-x(:,1);
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end

function [ obj ] =WFG4( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
A=30;
B=10;
C=0.35;

for i=1:c
    t1(:,i)=(1+cos((4*A+2)*pi*(0.5-(abs(x(:,i)-C))./(2*(floor(C-x(:,i))+C))))+4*B*(((abs(x(:,i)-C))./(2*(floor(C-x(:,i))+C))).^2))/(B+2);
end
for i=1:m-1
    y=t1(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    t2(:,i)=sum(y,2)/size(y,2);
end
y=t1(:,k+1:end);
t2(:,m)=sum(y,2)/size(y,2);
x=[];
for i=1:1
    x(:,i)=max(t2(:,m),ones(n,1)).*(t2(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t2(:,m),ones(n,1)).*(t2(:,i)-0.5)+0.5;
end
x(:,m)=t2(:,m);

h(:,1)=prod(sin(0.5*pi*x(:,1:m-1)),2);
for i=2:m-1
    h(:,i)=prod(sin(0.5*pi*x(:,1:m-i)),2).*(cos(0.5*pi*x(:,m-i+1)));
end
h(:,m)=cos(0.5*pi*x(:,1));
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end

function [ obj ] =WFG5( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
A=0.35;
B=0.001;
C=0.05;

for i=1:c
    t1(:,i)=1+(abs(x(:,i)-A)-B).*(floor(x(:,i)-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-x(:,i))*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end
for i=1:m-1
    y=t1(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    t2(:,i)=sum(y,2)/size(y,2);
end
y=t1(:,k+1:end);
t2(:,m)=sum(y,2)/size(y,2);
x=[];
for i=1:1
    x(:,i)=max(t2(:,m),ones(n,1)).*(t2(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t2(:,m),ones(n,1)).*(t2(:,i)-0.5)+0.5;
end
x(:,m)=t2(:,m);

h(:,1)=prod(sin(0.5*pi*x(:,1:m-1)),2);
for i=2:m-1
    h(:,i)=prod(sin(0.5*pi*x(:,1:m-i)),2).*(cos(0.5*pi*x(:,m-i+1)));
end
h(:,m)=cos(0.5*pi*x(:,1));
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end

function [ obj ] =WFG6( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
t1(:,1:k)=x(:,1:k);
t1(:,k+1:c)=abs(x(:,k+1:c)-0.35)./abs(floor(0.35-x(:,k+1:c))+0.35);

for i=1:m-1
    y=t1(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    A=k/(m-1);
    for j=1:size(y,2)
        s(:,j)=y(:,j);
        for p=0:A-2
            s(:,j)=s(:,j)+abs(y(:,j)-y(:,1+mod(j+p,size(y,2))));
        end
    end
    t2(:,i)=sum(s,2)/(size(y,2)/A*ceil(A/2)*(1+2*A-2*ceil(A/2)));
end
y=t1(:,k+1:end);
A=c-k;
for j=1:size(y,2)
    s(:,j)=y(:,j);
    for p=0:A-2
        s(:,j)=s(:,j)+abs(y(:,j)-y(:,1+mod(j+p,size(y,2))));
    end
end
t2(:,m)=sum(s,2)/(size(y,2)/A*ceil(A/2)*(1+2*A-2*ceil(A/2)));
x=[];
for i=1:1
    x(:,i)=max(t2(:,m),ones(n,1)).*(t2(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t2(:,m),ones(n,1)).*(t2(:,i)-0.5)+0.5;
end
x(:,m)=t2(:,m);

h(:,1)=prod(sin(0.5*pi*x(:,1:m-1)),2);
for i=2:m-1
    h(:,i)=prod(sin(0.5*pi*x(:,1:m-i)),2).*(cos(0.5*pi*x(:,m-i+1)));
end
h(:,m)=cos(0.5*pi*x(:,1));
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end

function [ obj ] =WFG7( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
A=0.98/49.98;
B=0.02;
C=50;
for i=1:k
    y=x(:,i+1:c);
    Y=sum(y,2)/size(y,2);
    for j=1:size(Y,1)
        if Y(j)<0.5
            u(j,1)=Y(j)*((C-B)*A)*2+B;
        else
            u(j,1)=(Y(j)-0.5)*(C-((C-B)*A+B))*2+(C-B)*A+B;
        end
    end
    v=A-(1-2*u).*abs(floor(0.5-u)+A);
    t1(:,i)=x(:,i).^(B+(C-B)*v);
end
t1(:,k+1:c)=x(:,k+1:c);

t2(:,1:k)=t1(:,1:k);
t2(:,k+1:c)=abs(t1(:,k+1:c)-0.35)./abs(floor(0.35-t1(:,k+1:c))+0.35);
for i=1:m-1
    y=t2(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    t3(:,i)=sum(y,2)/size(y,2);
end
y=t2(:,k+1:end);
t3(:,m)=sum(y,2)/size(y,2);
x=[];
for i=1:1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
x(:,m)=t3(:,m);

h(:,1)=prod(sin(0.5*pi*x(:,1:m-1)),2);
for i=2:m-1
    h(:,i)=prod(sin(0.5*pi*x(:,1:m-i)),2).*(cos(0.5*pi*x(:,m-i+1)));
end
h(:,m)=cos(0.5*pi*x(:,1));
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end

function [ obj ] =WFG8( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
A=0.98/49.98;
B=0.02;
C=50;
t1(:,1:k)=x(:,1:k);
for i=k+1:c
    y=x(:,1:i-1);
    Y=sum(y,2)/size(y,2);
    for j=1:size(Y,1)
        if Y(j)<0.5
            u(j,1)=Y(j)*((C-B)*A)*2+B;
        else
            u(j,1)=(Y(j)-0.5)*(C-((C-B)*A+B))*2+(C-B)*A+B;
        end
    end
    v=A-(1-2*u).*abs(floor(0.5-u)+A);
    t1(:,i)=x(:,i).^(B+(C-B)*v);
end

t2(:,1:k)=t1(:,1:k);
t2(:,k+1:c)=abs(t1(:,k+1:c)-0.35)./abs(floor(0.35-t1(:,k+1:c))+0.35);
for i=1:m-1
    y=t2(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    t3(:,i)=sum(y,2)/size(y,2);
end
y=t2(:,k+1:end);
t3(:,m)=sum(y,2)/size(y,2);
x=[];
for i=1:1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
x(:,m)=t3(:,m);

h(:,1)=prod(sin(0.5*pi*x(:,1:m-1)),2);
for i=2:m-1
    h(:,i)=prod(sin(0.5*pi*x(:,1:m-i)),2).*(cos(0.5*pi*x(:,m-i+1)));
end
h(:,m)=cos(0.5*pi*x(:,1));
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end

function [ obj ] =WFG9( m,k,z )
c=size(z,2);
n=size(z,1);
bu=2:2:2*c;
x=z./(ones(n,1)*bu);
A=0.98/49.98;
B=0.02;
C=50;
for i=1:c-1
    y=x(:,i+1:c);
    Y=sum(y,2)/size(y,2);
    for j=1:size(Y,1)
        if Y(j)<0.5
            u(j,1)=Y(j)*((C-B)*A)*2+B;
        else
            u(j,1)=(Y(j)-0.5)*(C-((C-B)*A+B))*2+(C-B)*A+B;
        end
    end
    v=A-(1-2*u).*abs(floor(0.5-u)+A);
    t1(:,i)=x(:,i).^(B+(C-B)*v);
end
t1(:,c)=x(:,c);

A=0.35;
B=0.001;
C=0.05;
for i=1:k
    t2(:,i)=1+(abs(t1(:,i)-A)-B).*(floor(t1(:,i)-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-t1(:,i))*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end
A=30;
B=95;
C=0.35;
for i=k+1:c
    t2(:,i)=(1+cos((4*A+2)*pi*(0.5-(abs(t1(:,i)-C))./(2*(floor(C-t1(:,i))+C))))+4*B*(((abs(t1(:,i)-C))./(2*(floor(C-t1(:,i))+C))).^2))/(B+2);
end

for i=1:m-1
    y=t2(:,(i-1)*k/(m-1)+1:i*k/(m-1));
    A=k/(m-1);
    for j=1:size(y,2)
        s(:,j)=y(:,j);
        for p=0:A-2
            s(:,j)=s(:,j)+abs(y(:,j)-y(:,1+mod(j+p,size(y,2))));
        end
    end
    t3(:,i)=sum(s,2)/(size(y,2)/A*ceil(A/2)*(1+2*A-2*ceil(A/2)));
end
y=t2(:,k+1:end);
A=c-k;
for j=1:size(y,2)
    s(:,j)=y(:,j);
    for p=0:A-2
        s(:,j)=s(:,j)+abs(y(:,j)-y(:,1+mod(j+p,size(y,2))));
    end
end
t3(:,m)=sum(s,2)/(size(y,2)/A*ceil(A/2)*(1+2*A-2*ceil(A/2)));
x=[];
for i=1:1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
for i=2:m-1
    x(:,i)=max(t3(:,m),ones(n,1)).*(t3(:,i)-0.5)+0.5;
end
x(:,m)=t3(:,m);

h(:,1)=prod(sin(0.5*pi*x(:,1:m-1)),2);
for i=2:m-1
    h(:,i)=prod(sin(0.5*pi*x(:,1:m-i)),2).*(cos(0.5*pi*x(:,m-i+1)));
end
h(:,m)=cos(0.5*pi*x(:,1));
s=2:2:2*m;
s=ones(n,1)*s;
obj=x(:,m)*ones(1,m)+s.*h;
end


