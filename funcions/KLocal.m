function [KL] = KLocal(L,E,A,v,ndofpn)
%v: vector barra sense normalitzar

v=v/norm(v); %normalitzat
Ke=zeros(ndofpn*2,ndofpn*2);
Ke(1,1)=1;
Ke(ndofpn+1,1)=-1;
Ke(1,ndofpn+1)=-1;
Ke(ndofpn+1,ndofpn+1)=1;

Ke=(E*A/L)*Ke;
x=[1; 0; 0];

u =cross(v,x); %Producte vectorial per saber eix de rotació

% if norm(cross(x,v))~=0
%     u=cross(x,v)./norm(cross(x,v));
% else
%     u=[0; 0; 1];   
% end
%theta = atan2(norm(cross(x,v)),dot(x,v)); %REVISAR!!
    %theta=theta(1);
    
theta=2*atan(norm(x*L-norm(x)*v)/norm(x*L+norm(x)*v));

%theta = acos(dot(x,v)/L); %Angle entre la barra i l'eix x

%theta=sign(sin(theta))*acos(cos(theta));

r=[cos(theta)+u(1)^2*(1-cos(theta))             u(1)*u(2)*(1-cos(theta))-u(3)*sin(theta)    u(1)*u(3)*(1-cos(theta))+u(2)*sin(theta)
   u(1)*u(2)*(1-cos(theta))+u(3)*sin(theta)     cos(theta)+u(2)^2*(1-cos(theta))            u(2)*u(3)*(1-cos(theta))-u(1)*sin(theta)
   u(1)*u(3)*(1-cos(theta))-u(2)*sin(theta)     u(2)*u(3)*(1-cos(theta))+u(1)*sin(theta)    cos(theta)+u(3)^2*(1-cos(theta))];

R = zeros(6,6);
R(1:3,1:3) = r;
R(4:6,4:6) = r;

KL = R\Ke*R;

end

