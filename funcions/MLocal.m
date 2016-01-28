function [ML] = MLocal(L,rho,A,v,alpha)

Mc=1/3*[2 0 0 1 0 0
        0 2 0 0 1 0
        0 0 2 0 0 1
        1 0 0 2 0 0
        0 1 0 0 2 0
        0 0 1 0 0 2];

Ml=eye(6);

M=1/2*rho*A*L*((1-alpha)*Mc+(alpha)*Ml);

x=[1; 0; 0];
v=v/norm(v);
u=cross(v,x); %Producte vectorial per saber eix de rotació

% if norm(cross(x,v))~=0
%     u=cross(x,v)./norm(cross(x,v));
% else
%     u=[0; 0; 1];   
% end

%     theta = atan2(dot(x',v)/norm(v)); %REVISAR!!
%     theta=theta(1);

%theta = acos(dot(x,v)/L); %Angle entre la barra i l'eix x

theta=2*atan(norm(x*L-norm(x)*v)/norm(x*L+norm(x)*v));

r=[cos(theta)+u(1)^2*(1-cos(theta))             u(1)*u(2)*(1-cos(theta))-u(3)*sin(theta)    u(1)*u(3)*(1-cos(theta))+u(2)*sin(theta)
   u(1)*u(2)*(1-cos(theta))+u(3)*sin(theta)     cos(theta)+u(2)^2*(1-cos(theta))            u(2)*u(3)*(1-cos(theta))-u(1)*sin(theta)
   u(1)*u(3)*(1-cos(theta))-u(2)*sin(theta)     u(2)*u(3)*(1-cos(theta))+u(1)*sin(theta)    cos(theta)+u(3)^2*(1-cos(theta))];

R = zeros(6,6);
R(1:3,1:3) = r;
R(4:6,4:6) = r;

ML=R\M*R;

end

