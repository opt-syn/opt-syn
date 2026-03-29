function s=syzfb(p, tol);
%function s=syzfb(p);
%
%Compute optimal convergence rate by bisection in algorithm synthesis. 
%
%The theory is exposed in C.W. Scherer, Ch. Ebenbauer, T. Holicki, 
%Optimization Algorithm Synthesis based on Integral Quadratic Constraints: A Tutorial, 
%62nd IEEE Conference on Decision and Control. 
%
%An extended version is available on arXiv under https://doi.org/10.48550/arXiv.2306.00565
%All references in the code are related to these paper.
%
%Bisection to optimize convergence rate rho.
%See syzf for required input arguments apart from p.rho.


if nargin < 2
    tol = 1e-4;
end

%smallest possible value of rho to satisfy conditions in Lemma 6
rhomin=max(abs(roots(p.alpha)));
if isempty(rhomin);rhomin=0;end;
rhoi=[rhomin,1];

ga=0;
s=[];
%bisection up to absolute accuracy of 1e-4
while rhoi(2)-rhoi(1)>tol;
    %select midpoint
    p.rho=(rhoi(2)+rhoi(1))/2;
    so=syzf(p);
    if so.ga<0 & norm(so.cd)>1e-3;
        %if successful, proceed with left interval 
        %keep s, and rho-values: interval in s.rhoi, last value in s.rho, minimal rho in s.rhomin
        rhoi(2)=p.rho;
        s=so;
        s.rhoi=rhoi;
        s.rho=p.rho;        
        s.rhomin=rhomin;
    else
        %if not successful, proceed with right interval
        rhoi(1)=p.rho;
    end;    
end;
   




