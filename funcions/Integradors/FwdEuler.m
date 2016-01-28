function [ y,timeComp] = FwdEuler( fun, y0, h, t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tic;
% nmax=tsim/h;
y(:,1)=y0;
for n=1:length(t)
%     t(n)=h*(n-1)+t0;
    y(:,n+1)=y(:,n)+h*fun(t(n),y(:,n));
end
timeComp=toc;
end

