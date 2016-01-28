clc;
clear all;
close all;
% Resolució P7 dinàmic

%Posició dels nodes i Matriu T de connectivitats

n=7; %número de costelles >=2
b=5; %envergadura
c=1.5; %corda
e=0.2; %espesor

[T, xpoints, npc]=NACA(n,b,c,e);

np=length(xpoints);

% Mòdul de young i àrea de les barres

E=6.5e10;
A=1e-3;
rho=2700;

%         node iDOF  Valor
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
vin(end-12,3)=0.5;
vin(end-15,3)=0.5;

 
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

% modes=zeros(ndofpn*nnode,length(vaps));
% modes(fFix,:)=veps;
modes=veps;
t=linspace(0,0.5,1000);
for iTemps=1:length(t)
    for iFreq=1:length(vaps)
%     t=linspace(0,2*pi/freq(iFreq),100);
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
subplot(2,1,2)
plot(t,Epot(:,:))
title('Energia potencial vs temps')
legend(legendCell)


figure(3)
b=bar(freq,mean(Etot'), 'basevalue', 10^-100);
set(gca,'YScale','log')
title('Energia total per freqüència')

% figure(4)
% loglog(freq,Etot(:,1),'kx')
% xlabel('Frequencia')
% ylabel('Energia')
  
% %Visualització amb GID
xnodes=xpoints;
conectivities=T;
npnod=nnode; %número total de nodes
nnpel=2; %número de nodes per element
ngaus=1;

% cd GID

% for iFreq=1:length(vaps)
%     mStr=int2str(iFreq);
%     folder=strcat('Mode',mStr);
%     mkdir(folder);
%     file_name=strcat(folder,'\P8_D');
%     t=linspace(0,2*pi/freq(iFreq),100);
%     modes=zeros(ndofpn*nnode,length(vaps));
%     modes(fFix,:)=veps;
%     mode(:,:,iFreq)=[reshape(modes(:,iFreq),3,[])]';
% 
%     for iTemps=1:length(t)
%         istep=iTemps;
% 
%         iter=istep;
% 
%         disp=mode(:,:,iFreq)*sin(freq(iFreq)*t(iTemps));
% 
%         %Post procés
%         ToGID (file_name,istep,xnodes,conectivities,nnpel,nelem,npnod);
%         ToGiD_post(file_name,iter,nnpel,ngaus,disp);
%     end
% end



%% Plot GiD
cd GID
file_name='P8_NACA';

%fix=[1 2 3 4 5 6];
fix=zeros(1,length(fixnodes)/ndofpn); %vector de 0
for i=1:length(fix)
        fix(i)=fixnodes(3*i-1,1);
end


dTotal=permute(sum(permute(x,[2 1 3])),[2 3 1]);
dFix=repmat(fixnodes(:,3),1,length(t));
dTemps=zeros(nnode*3, length(t));
dTemps(uFix,:)=dFix;
dTemps(fFix,:)=dTotal;


 for iTemps=1:length(t)
        disp=reshape(dTemps(:,iTemps),3,[])';
        %Post procés
        istep=iTemps;
        iter=istep;
        ToGID (file_name,istep,xnodes,conectivities,nnpel,nelem,npnod);
        ToGiD_post(file_name,iter,nnpel,ngaus,disp);
 end

% %% Integradors
% 
% I=eye(length(MNN));
% % dydt1=@(y2) I*y2; %x punt 1
% % dydt2=@(y1) -inv(MG)*KG*y1+inv(MG)*Fext; %x punt 2
% O=zeros(length(MNN));
% A=[O I; -inv(MNN)*KNN O];
% B=[zeros(length(MNN),1); inv(MNN)*Fext(:,3)];
% dydt=@(t,y) A*y+B;
% y0=[xin(:,3); vin(:,3)];
% 
% % t0=0;
% % tsim=0.1;
% h=t(end)/length(t); %2e-5
% 
% %Forward Euler
% [y1,timeComp1]=FwdEuler(dydt,y0,h,t);
% [Etot1,~,~]=Energia(y1,t,MNN,KNN);
% 
% %Backward Euler
% [y2,timeComp2]=BckEuler(dydt,y0,h,t);  
% [Etot2,~,~]=Energia(y2,t,MNN,KNN);
% 
% %RK4
% [y3,timeComp3]=RK4(dydt,y0,h,t); 
% [Etot3,~,~]=Energia(y3,t,MNN,KNN);
% 
% %ODE45
% [y4,t,timeComp4 ]=ODE45(dydt,y0,t );
% % tic;
% % [t4, y4,]=ode45(dydt,[t0 tsim],y0);
% % timeComp4=toc;
% [Etot4,~,~]=Energia(y4',t',MNN,KNN);
% 
% % %Mid-Point IRENE
% [y5,timeComp5]=MidPoint(dydt,y0,h,t);
% [Etot5,~,~]=Energia(y5,t,MNN,KNN);
% 
% for i=1:length(t)
%     Error1(i)=abs(Etot1(1)-Etot1(i))/Etot1(1);
%     Error2(i)=abs(Etot2(1)-Etot2(i))/Etot2(1);
%     Error3(i)=abs(Etot3(1)-Etot3(i))/Etot3(1);
%     Error4(i)=abs(Etot4(1)-Etot4(i))/Etot4(1);
%     Error5(i)=abs(Etot5(1)-Etot5(i))/Etot5(1);
% end
% 
% %Plots integradors
% figure(5)
% plot(t,Etot1);
% title('Energia total pel mètode Forward Euler')
% xlabel('Temps')
% ylabel('Energia')
% set(gca,'YScale','log')
% 
% 
% figure(6)
% plot(t,Error1);
% title('Error percentual en l´energia total pel mètode Forward Euler')
% xlabel('Temps')
% ylabel('Error')
% set(gca,'YScale','log')
% 
% figure(7)
% plot(t,Etot2);
% title('Energia total pel mètode Backward Euler')
% xlabel('Temps')
% ylabel('Energia')
% 
% figure(8)
% plot(t,Error2);
% title('Error percentual en l´energia total pel mètode Backward Euler')
% xlabel('Temps')
% ylabel('Error')
% 
% figure(9)
% plot(t,Etot3);
% title('Energia total pel mètode RK4')
% xlabel('Temps')
% ylabel('Energia')
% set(gca,'YScale','log')
% 
% figure(10)
% plot(t,Error3);
% title('Error percentual en l´energia total pel mètode RK4')
% xlabel('Temps')
% ylabel('Error')
% set(gca,'YScale','log')
% 
% figure(11)
% plot(t,Etot4);
% title('Energia total pel mètode Ode45 Matlab')
% xlabel('Temps')
% ylabel('Energia')
% 
% figure(12)
% plot(t,Error4);
% title('Error percentual en l´energia total pel mètode Ode45 Matlab')
% xlabel('Temps')
% ylabel('Error')
% 
% 
% figure(13)
% plot(t,Etot5);
% title('Energia total pel mètode Mid-Point IRENE')
% xlabel('Temps [s]')
% ylabel('Energia')
% 
% figure(14)
% plot(t,Error5);
% title('Error percentual en l´energia total pel mètode Mid-Point IRENE')
% xlabel('Temps[s]')
% ylabel('Error')
% 
% figure(15)
% subplot(2,3,1)
% plot(t,Error1)
% set(gca,'YScale','log')
% xlabel('Temps[s]')
% ylabel('Error Forward Euler')
% 
% subplot(2,3,2)
% plot(t,Error2)
% xlabel('Temps[s]')
% ylabel('Error Backward Euler')
% 
% subplot(2,3,3)
% plot(t,Error3)
% xlabel('Temps[s]')
% ylabel('Error RK4')
% set(gca,'YScale','log')
% 
% subplot(2,3,4)
% plot(t,Error4)
% xlabel('Temps[s]')
% ylabel('Error Ode45')
% 
% subplot(2,3,5)
% plot(t,Error5)
% xlabel('Temps[s]')
% ylabel('Error Mid Point')
% 
% mtit('Error dels diferents mètodes', 'fontsize',16,'xoff',-.0,'yoff',.025);
% 
%  
% time1=[timeComp1 timeComp2 timeComp3 timeComp4 timeComp5];
% m1=[1 2 3 4 5];
% figure
% bar(m1,time1);
% title('Temps computació de tots els mètodes')
% xlabel('Mètodes')
% ylabel('Temps de computació')
% 
% time2=[timeComp1 timeComp2 timeComp3 timeComp4];
% m2=[1 2 3 4];
% figure
% bar(m2,time2);
% title('Temps computació de tots els mètodes')
% xlabel('Mètodes')
% ylabel('Temps de computació')
