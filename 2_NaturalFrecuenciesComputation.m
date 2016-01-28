clc;
clear all;
close all;
% Resolució P7 dinàmic: modes de vibració 

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

E=6.5e10; %Pa
A=1e-4;
rho=2700; %kg/m^3

%         node iDOF  Valor
fixnodes = [1 1 0; 
            1 2 0;
            1 3 0;
            2 1 0;
            2 2 0;
            2 3 0;
            3 1 0;    
            3 2 0;
            3 3 0;
            4 1 0;
            4 2 0;
            4 3 0;
            ];

%    node iDOF  Valor
Fext = [5   1   0
        5   2   0
        5   3   0
        6   1   0
        6   2   0
        6   3   0
        7   1   0
        7   2   -1
        7   3   0
        8   1   0
        8   2   -1
        8   3   0];
    
ndofpn =3;
nnode=8;
KG=zeros(ndofpn*nnode,ndofpn*nnode);
MG=zeros(ndofpn*nnode,ndofpn*nnode);
nelem=22;

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
MNN=MG(fFix,fFix)
[v,D]=eig(MNN\KNN);
unitari=linspace(1,1,size(v,2));
vaps=D*unitari';
veps=v; %autovectors
freq=sqrt(abs(vaps)); %freqüències naturals (rad/s)
 
iFreq=1;

cd Dinamic
% mkdir(fprintf('Mode%1.0f',m))
% cd ..
% file_name=fprintf('Dinamic\Mode%1.0f\P7dinamic',m);

xnodes=xpoints;
conectivities=T;
npnod=nnode; %número total de nodes
nnpel=2; %número de nodes per element
ngaus=1;

for iFreq=1:length(vaps)
    mStr=int2str(iFreq);
    folder=strcat('Mode',mStr);
    mkdir(folder);
    file_name=strcat(folder,'\P7_D');
    t=linspace(0,2*pi/freq(iFreq),100);
    modes=zeros(ndofpn*nnode,length(vaps));
    modes(fFix,:)=veps;
    mode(:,:,iFreq)=[reshape(modes(:,iFreq),3,[])]';

    for iTemps=1:length(t)
        istep=iTemps;

        iter=istep;

        disp=mode(:,:,iFreq)*sin(freq(iFreq)*t(iTemps));

        %Post procés
        ToGID (file_name,istep,xnodes,conectivities,nnpel,nelem,npnod);
        ToGiD_post(file_name,iter,nnpel,ngaus,disp);
    end
end



