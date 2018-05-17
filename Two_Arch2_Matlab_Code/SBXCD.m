function [ NPOP ] = SBXCD( CA,DA,bu,bd,pc,n,m )
% Usage: [ NPOP ] = SBXCD( CA,DA,bu,bd,pc,n,m )
%
% Input:
% bu            -Upper Bound
% bd            -Lower Bound
% CA            -CA
% DA            -DA
% pc            -Crossover Probability
% n             -Population Scale
%
% Output: 
% NPOP          -Output Population by the Crossover between CA and DA
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
NPOP=[];
eta_c=15;
NC=size(CA,1);
ND=size(DA,1);
C=size(bu,2);
y=1;
for i=1:n
    r1=rand;
    if r1<=pc
        A=randperm(NC);
        k1=A(1);k2=A(2);
        x=dominated_relationship(CA(k1,:),CA(k2,:),m,C);
        if x==1
            k=k1;
        elseif x==2
            k=k2;
        else
            if rand>0.5
                k=k1;
            else
                k=k2;
            end
        end
        A=randperm(ND);
        y=A(1);
%         if k~=y %& d>0.001
            for j=1:C
                par1=DA(y,j);par2=CA(k,j);
                yd=bd(j);yu=bu(j);
                    r2=rand;
                    if r2<=0.5
                        y1=min(par1,par2);y2=max(par1,par2);
                        if (y1-yd)>(yu-y2)
                            beta=1+2*(yu-y2)/(y2-y1);
                        else
                            beta=1+2*(y1-yd)/(y2-y1);
                        end
                        expp=eta_c+1;beta=1/beta;alpha=2.0-beta^(expp);
                        r3=rand;
                        if r3<=1/alpha
                            alpha=alpha*r3;expp=1/(eta_c+1.0);
                            betaq=alpha^(expp);
                        else
                            alpha=1/(2.0-alpha*r3);expp=1/(eta_c+1);
                            betaq=alpha^(expp);
                        end
                        chld1=0.5*((y1+y2)-betaq*(y2-y1));
                        chld2=0.5*((y1+y2)+betaq*(y2-y1));   
                        aa=max(chld1,yd);
                        bb=max(chld2,yd);
                        if rand>0.5
                            NPOP(2*i-1,j)=min(aa,yu);
                            NPOP(2*i,j)=min(bb,yu);
                        else
                            NPOP(2*i,j)=min(aa,yu);
                            NPOP(2*i-1,j)=min(bb,yu);
                        end
                    else
                        NPOP(2*i-1,j)=par1;
                        NPOP(2*i,j)=par2;
                    end
            end
%         end
    end
    
end
end

