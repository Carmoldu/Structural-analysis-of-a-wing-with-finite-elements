% Resolució P7 dinàmic


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

E=1e9;
A=1e-4;
rho=6.3482;


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
freq=sqrt(abs(vaps)); %freqüències naturals
 

% Integradors
I=eye(length(MNN));

O=zeros(length(MNN));
A=[O I; -inv(MNN)*KNN O];
B=[zeros(length(MNN),1); inv(MNN)*Fext(:,3)];
dydt=@(t,y) A*y+B;
y0=[xin(:,3); vin(:,3)];

t0=0;
tsim=0.05;    

divsmax=12000000;
divsmin=10000000;
iterat=2;

Etot1=zeros(iterat,divsmax);

    
for j=1:iterat
    fprintf('Iteració %i \n',j)
    divs(j)=(j-1)*(divsmax-divsmin)/(iterat-1)+divsmin;
    
    %crear el vector detemps per cada iteració
    for k=1:divs(j)
        t(j,k)=(k-1)*(tsim)/(divs(j)-1)+t0;
    end
    
    h(j)=tsim/divs(j);
     

    %Forward Euler
    fprintf('\tMètode Forward Euler...')
    [y1,timeComp1(j)]=FwdEuler(dydt,y0,h(j),t(j,:));
    [Etot1(j,1:size(t,2)),~,~]=Energia(y1,t(j,:),MNN,KNN);
    fprintf('\t\tComplet!\tElapsed time: %fs\n',timeComp1(j))
    
    %Backward Euler
    fprintf('\tMètode Backward Euler...')
    [y2,timeComp2(j)]=BckEuler(dydt,y0,h(j),t(j,:));  
    [Etot2(j,1:size(t,2)),~,~]=Energia(y2,t(j,:),MNN,KNN);
    fprintf('\tComplet!\tElapsed time: %fs\n',timeComp2(j))
    
    %RK4
    fprintf('\tMètode RK4...')
    [y3,timeComp3(j)]=RK4(dydt,y0,h(j),t(j,:)); 
    [Etot3(j,1:size(t,2)),~,~]=Energia(y3,t(j,:),MNN,KNN);
    fprintf('\t\t\t\tComplet!\tElapsed time: %fs\n',timeComp3(j))

    %ODE45
    fprintf('\tMètode ODE45...')
    [y4,~,timeComp4(j)]=ODE45(dydt,y0,t(j,:));
    % tic;
    % [t4, y4,]=ode45(dydt,[t0 tsim],y0);
    % timeComp4=toc;
    [Etot4(j,1:size(t,2)),~,~]=Energia(y4',t(j,:)',MNN,KNN);
    fprintf('\t\t\t\tComplet!\tElapsed time: %fs\n',timeComp4(j))

    %Mid-Point IRENE
    fprintf('\tMètode Mid-Point...')
    [y5,timeComp5(j)]=MidPoint(dydt,y0,h(j),t(j,:));
    [Etot5(j,1:size(t,2)),~,~]=Energia(y5,t(j,:),MNN,KNN);
    fprintf('\t\t\tComplet!\tElapsed time: %fs\n',timeComp5(j))
    
    %BDF-2 CARLES
    fprintf('\tMètode BDF-2...')
    [y6,timeComp6(j)]=BDF2(dydt,y0,h(j),t(j,:));
    [Etot6(j,1:size(t,2)),~,~]=Energia(y6,t(j,:),MNN,KNN);
    fprintf('\t\t\tComplet!\tElapsed time: %fs\n',timeComp6(j))

    %Stormer-Verlet PERE ANTONI
    fprintf('\tMètode Stormer-Verlet...')
    [ x7,timeComp7(j) ] = Verlet( MNN,KNN,Fext,xin(:,3),vin(:,3),h(j),t(j,:));
    [ Etot7(j,1:size(t,2)), ~,~ ] = EnergiaX(x7,h(j),t(j,:),MNN,KNN);
    fprintf('\tComplet!\tElapsed time: %fs\n',timeComp7(j))
    

     for i=1:length(t(j,:))
         Error1(j,i)=abs(Etot1(j,1)-Etot1(j,i))/Etot1(j,1);
         Error2(j,i)=abs(Etot2(1)-Etot2(j,i))/Etot2(j,1);
         Error3(j,i)=abs(Etot3(1)-Etot3(j,i))/Etot3(j,1);
         Error4(j,i)=abs(Etot4(1)-Etot4(j,i))/Etot4(j,1);
         Error5(j,i)=abs(Etot5(1)-Etot5(j,i))/Etot5(j,i);
         Error6(j,i)=abs(Etot6(1)-Etot6(j,i))/Etot6(j,i);
         Error7(j,i)=abs(Etot7(1)-Etot7(j,i))/Etot7(j,i);
     end
    fprintf('\t\t\t\t\t\t\t\t\t\t\tTotal Elapsed time: %fs\n\n',timeComp1(j)+timeComp2(j)+timeComp3(j)+timeComp4(j)+timeComp5(j)+timeComp6(j))
end

Error1(Error1==0) = nan ;
Error2(Error2==0) = nan ;
Error3(Error3==0) = nan ;
Error4(Error4==0) = nan ;
Error5(Error5==0) = nan ;
Error6(Error6==0) = nan ;
Error7(Error7==0) = nan ;

%Plots integradors
figure(6)
plot(t',Error1');
title('Error percentual en l´energia total pel mètode Forward Euler')
xlabel('Temps')
ylabel('Error')
set(gca,'YScale','log')
legendCell=strcat('h [s]= ',strtrim(cellstr(num2str(h'))));
legend(legendCell)
 
figure(8)
plot(t',Error2');
title('Error percentual en l´energia total pel mètode Backward Euler')
xlabel('Temps')
ylabel('Error')
legendCell=strcat('h [s]= ',strtrim(cellstr(num2str(h'))));
legend(legendCell)

figure(10)
plot(t',Error3');
title('Error percentual en l´energia total pel mètode RK4')
xlabel('Temps')
ylabel('Error')
set(gca,'YScale','log')
legendCell=strcat('h [s]= ',strtrim(cellstr(num2str(h'))));
legend(legendCell)

figure(12)
plot(t',Error4');
title('Error percentual en l´energia total pel mètode Ode45 Matlab')
xlabel('Temps')
ylabel('Error')
legendCell=strcat('h [s]= ',strtrim(cellstr(num2str(h'))));
legend(legendCell)

figure(14)
plot(t',Error5');
title('Error percentual en l´energia total pel mètode Mid-Point IRENE')
xlabel('Temps[s]')
ylabel('Error')
legendCell=strcat('h [s]= ',strtrim(cellstr(num2str(h'))));
legend(legendCell)

figure(13)
plot(t(10,:)',Error6');
title('Error percentual en l´energia total pel mètode BDF-2 CARLES')
xlabel('Temps[s]')
ylabel('Error')
legendCell=strcat('h [s]= ',strtrim(cellstr(num2str(h'))));
legend(legendCell)
set(gca,'YScale','log')

figure(11)
plot(t',Error7');
title('Error percentual en l´energia total pel mètode Stormer-Verlet PERE ANTONI')
xlabel('Temps[s]')
ylabel('Error')
legendCell=strcat('h [s]= ',strtrim(cellstr(num2str(h'))));
legend(legendCell)
set(gca,'YScale','log')


figure(15)
subplot(2,4,1)
plot(t',Error1')
set(gca,'YScale','log')
xlabel('Temps[s]')
ylabel('Error Forward Euler')

subplot(2,4,2)
plot(t',Error2')
xlabel('Temps[s]')
ylabel('Error Backward Euler')

subplot(2,4,3)
plot(t',Error3')
xlabel('Temps[s]')
ylabel('Error RK4')
set(gca,'YScale','log')

subplot(2,4,4)
plot(t',Error4')
xlabel('Temps[s]')
ylabel('Error Ode45')

subplot(2,4,5)
plot(t',Error5')
xlabel('Temps[s]')
ylabel('Error Mid Point')

subplot(2,4,6)
plot(t(10,:)',Error6')
xlabel('Temps[s]')
ylabel('Error BDF-2')
set(gca,'YScale','log')

subplot(2,4,7)
plot(t(10,:)',Error7')
xlabel('Temps[s]')
ylabel('Error Stormer-Verlet')
set(gca,'YScale','log')
mtit('Error dels diferents mètodes', 'fontsize',16,'xoff',-.0,'yoff',.025);

legendCell=strcat('h [s]= ',strtrim(cellstr(num2str(h'))));
legend(legendCell)

 
figure(16)
hold on
plot(h,timeComp1,'color','b')
plot(h,timeComp2,'color','m')
plot(h,timeComp3,'color','c')
plot(h,timeComp4,'color','r')
plot(h,timeComp5,'color','g')
plot(h,timeComp6,'color','b')
plot(h,timeComp7,'color','m')

legend('Forward Euler')
title('Temps de computació per a tots els mètodes explícits per a diferents h')
xlabel('h [s]')
ylabel('Temps computació [s]')
set(gca,'YScale','log')
