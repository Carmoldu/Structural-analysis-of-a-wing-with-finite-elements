function [ y, timeComp  ] = BDF2( fun, y0, h, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tic;
% nmax=tsim/h;
y(:,1)=y0;
alpha=3/2;

options=optimset('Display','off');
for n=1:length(t)-1
    if n==1
        faux=@(yplus) y(:,n)+h*fun(t(n)+h,yplus)-yplus; %funcio auxiliar
        y(:,n+1)=fsolve(faux,y(:,n),options);
    else 
        faux=@(yplus) h*fun(t(n)+h,yplus)-(1-alpha)*(y(:,n)-y(:,n-1))-alpha*(yplus-y(:,n)); %funcio auxiliar
        y(:,n+1)=fsolve(faux,y(:,n),options);
    end
end  
timeComp=toc;
end

