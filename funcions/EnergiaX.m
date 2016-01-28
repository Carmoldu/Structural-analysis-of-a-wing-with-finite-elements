function [  Etot, Ecin, Epot ] = EnergiaX(X,h,t,MNN,KNN)
    x=X;
    v=zeros(size(x));
    for i=1:length(t)-1
        v(:,i)=(x(:,i+1)-x(:,i))/h;
    end
    v(:,end)=v(:,end-1);
    for iTemps=1:length(t)      
            Ecin(iTemps)=1/2*v(:,iTemps)'*MNN* v(:,iTemps); %Energia cinètica
            Epot(iTemps)=1/2*x(:,iTemps)'*KNN*x(:,iTemps); %Energia potencial
            Etot(iTemps)=Ecin(iTemps)+Epot(iTemps);        
    end
end

