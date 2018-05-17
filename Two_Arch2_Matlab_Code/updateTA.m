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

function [ NCA,NDA] = updateTA(CPOP,CA,NPOP1,DA,c,m,n,p)
% Usage: [ NCA,NDA] = updateTA(CPOP,CA,NPOP1,DA,c,m,n,p)
%
% Input:
% CPOP          -Non-dominated Solution Offspring Set
% CA            -CA
% DA            -DA
% NPOP1         -Offspring Set
% m             -No. of Objectives
% c             -No. of Decision Variables
% n             -Population Scale
% p             -Lp-norm
%
% Output: 
% NCA           -Updated CA
% NDA           -Updated DA
%

%Update CA
nc=size(CPOP,1);
P=[CA;NPOP1];
P=find_unique( P,c,m);
F = Fitness( P,c,m );
[ NCA,F] = Environment_selection( P,c,m,F,n);
J=[];
for i=1:nc
    flag=0;
    np=0;
    I=[];
    for j=1:size(DA,1)
        x=dominated_relationship(CPOP(i,:),DA(j,:),m,c);
        if x==2
            flag=1;
            break;
        elseif x==1
            np=np+1;
            I=[I;j];
        elseif x==3
            flag=1;
            break;
        end
    end
    if size(I,1)~=0
        DA(I,:)=[];
    end
    if flag==0
        J=[J,i];
    end
end
%Update DA
DA=[DA;CPOP(J,:)];
if size(DA,1)>n
    NDA= selectDAlp( DA,n,c,m,p);
else
    NDA=DA;
end

end

function [ NPOP ] = find_unique( POP,c,m)
% Usage: [ NPOP ] = find_unique( POP,c,m)
%
% Input:
% POP           -Population
% m             -No. of Objectives
% c             -No. of Decision Variables
%
% Output: 
% NPOP          -Non-duplicated Population
%
n=size(POP,1);
i=1;
while i<=size(POP,1)
    flag=0;
    j=i+1;
    while j<=size(POP,1)
        x=dominated_relationship(POP(i,:),POP(j,:),m,c);
        if x==3
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

function [ F ] = Fitness( POP,c,m )
% Usage: [ F ] = Fitness( POP,c,m )
%
% Input:
% POP           -Population
% m             -No. of Objectives
% c             -No. of Decision Variables
%
% Output: 
% F             -Fitness of IBEA
n=size(POP,1);
CPOP=POP(:,c+1:c+m);
bu=max(CPOP);
bd=min(CPOP);
CPOP=(CPOP-repmat(bd,n,1))./(repmat(bu,n,1)-repmat(bd,n,1));
F=zeros(n,1);
for i=1:n
    t=CPOP(:,1:m);
    t(i,:)=[];
    t1=repmat(CPOP(i,1:m),n-1,1);
    t2=max(t-t1,[],2);
    c=max(abs(t2));
    F(i)=sum(0-exp(0-t2/(c*0.05)));
end
end

function [ NPOP,F ] = Environment_selection( POP,c,m,F,n)
% Usage: [ NPOP,F ] = Environment_selection( POP,c,m,F,n)
%
% Input:
% POP           -Population
% m             -No. of Objectives
% c             -No. of Decision Variables
% n             -Population Scale
% F             -Fitness of IBEA
%
% Output: 
% F             -Fitness of IBEA
% NPOP          -Selected Solution Set
CPOP=POP(:,c+1:c+m);
bu=max(CPOP);
bd=min(CPOP);
N=size(CPOP,1);
CPOP=(CPOP-repmat(bd,N,1))./(repmat(bu,N,1)-repmat(bd,N,1));

C=zeros(N,1);
for i=1:N
    t=CPOP(:,1:m);
    t(i,:)=[];
    t1=repmat(CPOP(i,1:m),N-1,1);
    t2=max(t-t1,[],2);
    C(i)=max(abs(t2));
end

while size(POP,1)>n
    [t,I]=min(F);
    t1=repmat(CPOP(I,1:m),size(POP,1)-1,1);
    CPOP(I,:)=[];
    t2=max(t1-CPOP,[],2);
    POP(I,:)=[];
    F(I)=[];
    F=F+exp(0-t2/(C(I)*0.05));
    C(I)=[];
end
NPOP=POP;
end


function [ NPOP] = selectDAlp( POP,n,c,m,p)
% Usage: [ NPOP] = selectDAlp( POP,n,c,m,p)
%
% Input:
% POP           -Overflowed DA
% m             -No. of Objectives
% c             -No. of Decision Variables
% n             -Population Scale
% p             -Lp-norm
%
% Output: 
% NPOP          -Updated DA
NPOP=[];
[A,I1]=max(POP(:,c+1:c+m));
[A,I2]=min(POP(:,c+1:c+m));
I=[I1,I2];
I=unique(I);
%Bunderay Points
NPOP=[NPOP;POP(I,:)];
POP(I,:)=[];
while size(NPOP,1)<n & size(POP,1)>0
    d=zeros(size(POP,1),1);
    for i=1:size(POP,1)
        dis=sum((abs(NPOP(:,c+1:c+m)-ones(size(NPOP,1),1)*POP(i,c+1:c+m))).^p,2).^(1/p);
        d(i)=min(dis);
    end
    [A,I]=max(d);
    Imin=find(d<0.001);
    NPOP=[NPOP;POP(I,:)];
    POP([I;Imin],:)=[];
end
end
