clc;
clear all;
close all;
% Resolució dinàmic

% Posició dels nodes i Matriu T de connectivitats

xpoints = [0               0               1
           0               0               0
           0               1               1
           0               1               0
           2            0.25            0.75
           2            0.25            0.25
           2            0.75            0.75
           2            0.75            0.25];

    % Matriu T de connectivitats

T = [1 5
     5 7
     7 3
     8 6
     6 5
     8 7
     4 8
     2 6
     7 6
     8 5
     1 7
     5 3
     6 4
     2 8
     7 4
     8 3
     6 1
     5 2
     3 1
     3 4
     4 2
     2 1];

% Mòdul de young i àrea de les barres

E=6.5e10;
A=1e-4;
rho=2700;


%         node iDOF  Valor
npc=4;
np=length(xpoints);

fixnodes=zeros((npc*3),3);

for i=1:npc
    for j=1:3
        k=3*(i-1)+j;
        fixnodes(k,1)=i;
        fixnodes(k,2)=j;
    end
end
%    node iDOF  Valor
Fext=zeros((np-npc)*3,3);
for i=1:np-npc
    for j=1:3
        k=3*(i-1)+j;
        Fext(k,1)=i+npc;
        Fext(k,2)=j;
    end
end

% Fext = [5   1   0
%         5   2   0
%         5   3   0
%         6   1   0
%         6   2   0
%         6   3   0
%         7   1   0
%         7   2   -1
%         7   3   0
%         8   1   0
%         8   2   -1
%         8   3   0];
    
% vel inicial
vin=zeros((np-npc)*3,3);
for i=1:np-npc
    for j=1:3
        k=3*(i-1)+j;
        vin(k,1)=i+npc;
        vin(k,2)=j;
    end
end
 % posicio inicial
xin=zeros((np-npc)*3,3);
for i=1:np-npc
    for j=1:3
        k=3*(i-1)+j;
        xin(k,1)=i+npc;
        xin(k,2)=j;
    end
end

%Definició dels punts inicials no-nuls
vin(end,3)=0.5;
vin(end-3,3)=0.5;
vin(end-6,3)=0.5;
vin(end-9,3)=0.5;

ndofpn =3;
nnode=np;
KG=zeros(ndofpn*nnode,ndofpn*nnode);
MG=zeros(ndofpn*nnode,ndofpn*nnode);
nelem=length(T);

alpha=0; %Per la matriu masses local

%B) Càlcul matriu global
for ielem=1:nelem
    v(:,ielem) = xpoints(T(ielem,2),:)-xpoints(T(ielem,1),:); %vector barra
    L(ielem)=norm(v(:,ielem)); %longitud barra
    KL(:,:,ielem)=KLocal(L(ielem),E,A,v(:,ielem),ndofpn);
    ML(:,:,ielem)=MLocal(L(ielem),rho,A,v(:,ielem),alpha);
    for inode=1:2 %numero de nodes per barra
        for jnode=1:2 %numero de nodes per barra
            for idofpn=1:ndofpn
                for jdofpn=1:ndofpn
                    iG=idofpn+(T(ielem,inode)-1)*ndofpn;
                    jG=jdofpn+(T(ielem,jnode)-1)*ndofpn;
                    iL=idofpn+ndofpn*(inode-1);
                    jL=jdofpn+ndofpn*(jnode-1);
                    KG(iG,jG)=KG(iG,jG)+KL(iL,jL,ielem);
                    MG(iG,jG)=MG(iG,jG)+ML(iL,jL,ielem);
                end
            end
        end
    end
end 

uFix=fixnodes(:,2)+ndofpn*(fixnodes(:,1)-1);

fFix=setdiff(1:ndofpn*nnode,uFix);
KNN=KG(fFix,fFix);
KDD=KG(uFix,uFix);
KDN=KG(uFix,fFix);
KND=KG(fFix,uFix);
fN=Fext(:,3);
uD=fixnodes(:,3);

uN=KNN\(fN-KND*uD);
fD=KDD*uD+KDN*uN;
f=zeros(ndofpn*nnode,1);
u=zeros(ndofpn*nnode,1);
f(fFix)=fN;
f(uFix)=fD; %
u(fFix)=uN;
u(uFix)=uD; %

%Càlcul freqüències naturals

MNN=MG(fFix,fFix);
[v,D]=eig(MNN\KNN);
unitari=linspace(1,1,size(v,2));
vaps=D*unitari';
veps=v; %autovectors
freq=sqrt(abs(vaps)); %freqüències naturals (rad/s)
freqHz=freq/(2*pi); %Freqüències naturals (Hz)
 
% Calcul a i b
a=veps\xin(:,3);
b=veps\vin(:,3)./freq;


modes=veps;
t=linspace(0,0.05,1e4); %1e4
for iTemps=1:length(t)
    for iFreq=1:length(vaps)
        x(:,iFreq,iTemps)=modes(:,iFreq)*(a(iFreq)*cos(freq(iFreq)*t(iTemps))+b(iFreq)*sin(freq(iFreq)*t(iTemps)));
        v(:,iFreq,iTemps)=modes(:,iFreq)*freq(iFreq)*(b(iFreq)*cos(freq(iFreq)*t(iTemps))-a(iFreq)*sin(freq(iFreq)*t(iTemps)));
        Ecin(iFreq,iTemps)=1/2*v(:,iFreq,iTemps)'*MNN* v(:,iFreq,iTemps); %Energia cinètica
        Epot(iFreq,iTemps)=1/2*x(:,iFreq,iTemps)'*KNN*x(:,iFreq,iTemps); %Energia potencial
        Etot(iFreq,iTemps)=Ecin(iFreq,iTemps)+Epot(iFreq,iTemps);
    end
end

%suma x i v per totes les freqüències (només depèn del temps
for i=1:length(t)
    for j=1:length(fFix)
        X(j,i)=sum(x(j,:,i));
        V(j,i)=sum(v(j,:,i));
    end
    EC.TOT(i)=sum(Ecin(:,i));
    EP.TOT(i)=sum(Epot(:,i));
    ET.TOT(i)=sum(Etot(:,i));
end
Y=[X;V]; %Vector de desplaçaments i velocitats concatenats per comparar amb els integradors

% PLOTs
figure(1)
semilogy(t,Etot(:,:))
title('Energia total vs temps')
xlabel('t (s)');
ylabel('Energia total');
%legendCell=strcat('Freq [MHz]= ',strtrim(cellstr(num2str(floor(vaps*10^(-6)*100)/100))));
legendCell=strcat('Freq [Hz]= ',strtrim(cellstr(num2str(floor(freqHz)))));
legend(legendCell)

figure(2)
subplot(2,1,1)
plot(t,Ecin(:,:))
title('Energia cinètica vs temps')
xlabel('t (s)');
ylabel('Energia cinètica');
subplot(2,1,2)
plot(t,Epot(:,:))
title('Energia potencial vs temps')
xlabel('t (s)');
ylabel('Energia potencial');
legend(legendCell)

figure(3)
b=bar(1:length(freq),mean(Etot'));%, 'basevalue', 10^-100);
%set(gca,'YScale','log')
title('Energia total per freqüència')
xlabel('Mode');
ylabel('Energia total');

figure(4)
loglog(freqHz,Etot(:,1),'kx')
xlabel('Frequencia (Hz)')
ylabel('Energia')

% figure
% plot(t,EC.TOT,'r',t,EP.TOT,'k',t,ET.TOT)
% title('Energies del sistema')
% xlabel('t (s)');
% ylabel('Energia');

xnodes=xpoints;
conectivities=T;
npnod=nnode; %número total de nodes
nnpel=2; %número de nodes per element
ngaus=1;

%% Plot GID
% cd GID
% file_name='P8';
% 
% %fix=[1 2 3 4 5 6];
% fix=zeros(1,length(fixnodes)/ndofpn); %vector de 0
% for i=1:length(fix)
%         fix(i)=fixnodes(3*i-1,1);
% end
% 
% dTotal=permute(sum(permute(x,[2 1 3])),[2 3 1]);
% dFix=repmat(fixnodes(:,3),1,length(t));
% dTemps=zeros(nnode*3, length(t));
% dTemps(uFix,:)=dFix;
% dTemps(fFix,:)=dTotal;
% 
%  for iTemps=1:length(t)
%         disp=reshape(dTemps(:,iTemps),3,[])';
%         %Post procés
%         istep=iTemps;
%         iter=istep;
%         ToGID (file_name,istep,xnodes,conectivities,nnpel,nelem,npnod);
%         ToGiD_post(file_name,iter,nnpel,ngaus,disp);
%  end
%% Integradors

I=eye(length(MNN));
% dydt1=@(y2) I*y2; %x punt 1
% dydt2=@(y1) -inv(MG)*KG*y1+inv(MG)*Fext; %x punt 2
O=zeros(length(MNN));
A=[O I; -inv(MNN)*KNN O];
B=[zeros(length(MNN),1); inv(MNN)*Fext(:,3)];
dydt=@(t,y) A*y+B;
y0=[xin(:,3); vin(:,3)];

% t0=0;
% tsim=0.1;
h=t(end)/length(t); 

%Forward Euler
[y1,timeComp1]=FwdEuler(dydt,y0,h,t);
[Etot1,~,~]=Energia(y1,t,MNN,KNN);

%Backward Euler
[y2,timeComp2]=BckEuler(dydt,y0,h,t);  
[Etot2,~,~]=Energia(y2,t,MNN,KNN);

%RK4
[y3,timeComp3]=RK4(dydt,y0,h,t); 
[Etot3,~,~]=Energia(y3,t,MNN,KNN);

%ODE45
[y4,t4,timeComp4 ]=ODE45(dydt,y0,t );
y4=y4';
t4=t4';
[Etot4,~,~]=Energia(y4,t,MNN,KNN);

%Mid-Point IRENE
[y5,timeComp5]=MidPoint(dydt,y0,h,t);
[Etot5,~,~]=Energia(y5,t,MNN,KNN);

%BDF-2 CARLES
[y6,timeComp6]=BDF2(dydt,y0,h,t);
[Etot6,~,~]=Energia(y6,t,MNN,KNN);

%Stormer-Verlet PERE ANTONI
[ x7,timeComp7 ] = Verlet( MNN,KNN,Fext,xin(:,3),vin(:,3),h,t);
[ Etot7, Ecin7, Epot7 ] = EnergiaX(x7,h,t,MNN,KNN);

for i=1:length(t)
    %Error energia
    Error1(i)=abs(Etot1(1)-Etot1(i))/Etot1(1);
    Error2(i)=abs(Etot2(1)-Etot2(i))/Etot2(1);
    Error3(i)=abs(Etot3(1)-Etot3(i))/Etot3(1);
    Error4(i)=abs(Etot4(1)-Etot4(i))/Etot4(1);
    Error5(i)=abs(Etot5(1)-Etot5(i))/Etot5(1);
    Error6(i)=abs(Etot6(1)-Etot6(i))/Etot6(i);
    Error7(i)=abs(Etot7(1)-Etot7(i))/Etot7(i);
    
    %Error y analitic vs y integradors
    e1(i)=norm(Y(:,i)-y1(:,i))/norm(Y(:,i));
    e2(i)=norm(Y(:,i)-y2(:,i))/norm(Y(:,i));
    e3(i)=norm(Y(:,i)-y3(:,i))/norm(Y(:,i));
    e4(i)=norm(Y(:,i)-y4(:,i))/norm(Y(:,i));
    e5(i)=norm(Y(:,i)-y5(:,i))/norm(Y(:,i));
    e6(i)=norm(Y(:,i)-y6(:,i))/norm(Y(:,i));
    e7(i)=norm(Y(1:end/2,i)-x7(:,i))/norm(Y(1:end/2,i)); %SV es basa en la posició
    
end

%Plots integradors
figure
plot(t,Etot,t,Etot1,'r',t,Etot2,'c',t,Etot3,'k',t,Etot4,'g',t,Etot5,'m',t,Etot6,'y',t,Etot7,'b');
legend('Analític','FwdEuler','BckEuler','RK4','Ode45','Midpoint','BDF-2','Verlet')
title('Energies totals vs temps')
xlabel('Temps')
ylabel('Energia')

figure
semilogy(t,Error1,'r',t,Error2,'c',t,Error3,'k',t,Error4,'g',t,Error5,'m',t,Error6,'y',t,Error7);
legend('FwdEuler','BckEuler','RK4','Ode45','Midpoint','BDF-2','Verlet')
title('Error energia vs temps')
xlabel('Temps')
ylabel('Error')

figure
semilogy(t,e1,'r',t,e2,'c',t,e3,'k',t,e4,'g',t,e5,'m',t,e6,'y',t,e7);
legend('FwdEuler','BckEuler','RK4','Ode45','Midpoint','BDF-2','Verlet')
title('Error analític/integrador vs temps')
xlabel('Temps')
ylabel('Error')

% figure
% plot(t,Etot1);
% title('Forward Euler')
% xlabel('Temps')
% ylabel('Energia')
% 
% figure
% semilogy(t,Error1);
% title('Forward Euler')
% xlabel('Temps')
% ylabel('Error')
% 
% figure
% plot(t,Etot2);
% title('Backward Euler')
% xlabel('Temps')
% ylabel('Energia')
% 
% figure
% semilogy(t,Error2);
% title('Backward Euler')
% xlabel('Temps')
% ylabel('Error')
% 
% figure
% plot(t,Etot3);
% title('RK4')
% xlabel('Temps')
% ylabel('Energia')
% 
% figure
% semilogy(t,Error3);
% title('RK4')
% xlabel('Temps')
% ylabel('Error')
% 
% figure
% plot(t,Etot4);
% title('Ode45 Matlab')
% xlabel('Temps')
% ylabel('Energia')
% 
% figure
% semilogy(t,Error4);
% title('Ode45 Matlab')
% xlabel('Temps')
% ylabel('Error')
% 
% figure
% plot(t,Etot5);
% title('Mid-Point IRENE')
% xlabel('Temps')
% ylabel('Energia')
% 
% figure
% semilogy(t,Error5);
% title('Mid-Point IRENE')
% xlabel('Temps')
% ylabel('Error')
 
time1=[timeComp1 timeComp2 timeComp3 timeComp4 timeComp5 timeComp6 timeComp7];
m1=[1 2 3 4 5 6 7];
figure
bar(m1,time1);
title('Temps computació de tots els mètodes')
xlabel('Mètodes')
ylabel('Temps de computació')

