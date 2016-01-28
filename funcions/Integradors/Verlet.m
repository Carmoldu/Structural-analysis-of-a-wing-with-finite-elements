function [ x,tcomp ] = Verlet( MNN,KNN,Fext,x0,v0,h,t)
    Av=@(t,x) MNN\(Fext(:,3)-KNN*x);
    tic
    x(:,1)=x0;
    x(:,2)=x0+v0*h+0.5*Av(t(1),x0)*h^2;
    for i=3:length(t)
        x(:,i)=2*x(:,i-1)-x(:,i-2)+Av(t(i-1),x(:,i-1))*h^2;
    end
    tcomp=toc;

end

