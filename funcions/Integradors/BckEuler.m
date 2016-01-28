function [ y, timeComp  ] = BckEuler( fun, y0, h, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tic;
% nmax=tsim/h;
y(:,1)=y0;
options=optimset('Display','off');
for n=1:length(t)-1
%     t(n)=h*(n-1)+t0;
    faux=@(yplus) y(:,n)+h*fun(t(n)+h,yplus)-yplus; %funcio auxiliar
    y(:,n+1)=fsolve(faux,y(:,n),options);
end  
timeComp=toc;
end

