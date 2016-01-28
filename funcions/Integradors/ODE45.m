function [ y, T, timeComp  ] = ODE45( fun, y0, t )

tic;
[T, y]=ode45(fun,t,y0); %t: time
timeComp=toc;

end

