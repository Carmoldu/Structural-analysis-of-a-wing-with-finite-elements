clear;
clc;
close all;

%Problema 7 Estàtic
    %A) DADES

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
nelem=22;

%B) Càlcul matriu global
for ielem=1:nelem
    v(:,ielem) = xpoints(T(ielem,2),:)-xpoints(T(ielem,1),:); %vector barra
    L(ielem)=norm(v(:,ielem)); %longitud barra
    %w(:,ielem) = v(:,ielem)./L(ielem); %vector barra normalitzat
    KL(:,:,ielem) = KLocal(L(ielem),E,A,v(:,ielem),ndofpn);
    for inode=1:2 %numero de nodes per barra
        for jnode=1:2 %numero de nodes per barra
            for idofpn=1:ndofpn
                for jdofpn=1:ndofpn
                    iG=idofpn+(T(ielem,inode)-1)*ndofpn;
                    jG=jdofpn+(T(ielem,jnode)-1)*ndofpn;
                    iL=idofpn+ndofpn*(inode-1);
                    jL=jdofpn+ndofpn*(jnode-1);
                    KG(iG,jG)=KG(iG,jG)+KL(iL,jL,ielem);
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
f(uFix)=fD
u(fFix)=uN;
u(uFix)=uD
    
istep=1;
file_name='Estatic\P7estatic';
xnodes=xpoints;
conectivities=T;
npnod=nnode; %número total de nodes
nnode=2; %número de nodes per element
iter=istep;
disp=u; %desplaçaments: [ux; uy; uz]
ngaus=1;
disp=[reshape(u,3,[])]';

%Post procés
ToGID (file_name,istep,xnodes,conectivities,nnode,nelem,npnod);
ToGiD_post(file_name,iter,nnode,ngaus,disp);
