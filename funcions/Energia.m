function [ Etot, Ecin, Epot] = Energia(Y,t,MNN,KNN)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x=Y(1:size(Y,1)/2,:);
v=Y(size(Y,1)/2+1:size(Y,1),:);
for iTemps=1:length(t)      
        Ecin(iTemps)=1/2*v(:,iTemps)'*MNN* v(:,iTemps); %Energia cinètica
        Epot(iTemps)=1/2*x(:,iTemps)'*KNN*x(:,iTemps); %Energia potencial
        Etot(iTemps)=Ecin(iTemps)+Epot(iTemps);        
end
end

