function [ y, timeComp  ] = RK4( fun, y0, h, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tic;
% nmax=tsim/h;
y(:,1)=y0;
for n=1:length(t)
%     t(n)=h*(n-1)+t0;
    K1=fun(t(n),y(:,n));
    K2=fun(t(n)+h/2,y(:,n)+h/2*K1);
    K3=fun(t(n)+h/2,y(:,n)+h/2*K2);
    K4=fun(t(n)+h,y(:,n)+h*K3);
    
    y(:,n+1)=y(n)+h/6*(K1+2*K2+2*K3+K4);
end  
timeComp=toc;


end
