function [ y,timeComp   ] =MidPoint( fun, y0, h, t )

tic;
% nmax=tsim/h;
y(:,1)=y0;
options=optimset('Display','off');
for n=1:length(t)
%     t(n)=h*(n-1)+t0;
    faux=@(yplus) y(:,n)+h*fun(t(n)+h/2,1/2*(yplus+y(:,n)))-yplus; %funcio auxiliar
    y(:,n+1)=fsolve(faux,y(:,n),options);
end  
timeComp=toc;

end

