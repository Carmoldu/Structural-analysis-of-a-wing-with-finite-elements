function [ T,xpoints,npc ] = NACA(n,b,c,e)
% n=2; %número de costelles >=2
% b=5; %envergadura
% c=1; %corda
% e=0.2; %espesor
npc=6; %número de nodes per costella
epc=(npc-4)*2+5; %número d'elements per costella
dim=3; %dimensions
dc=b/(n-1); %distància entre costelles
xpoints=zeros(npc*n,dim);
perfilNACA= [0 0
            c/4 -e/2
            c/4 e/2
            c/2 -e/2
            c/2 e/2
            c 0
            ];
ConectivitatsNACA=[1 2
                   1 3
                   2 3
                   2 4
                   3 5
                   2 5
                   4 5
                   4 6
                   5 6];
           
conectivitatsEstructura=[
                        1 1 
                        1 2
                        1 3
                        2 2
                        2 3
                        2 4
                        3 3
                        3 5
                        4 4
                        4 5
                        4 6
                        5 5
                        5 6
                        6 6
                        ];
            
for i=1:n
    xpoints(npc*(i-1)+1:npc*i,2)= (i-1)*dc ;
    xpoints(npc*(i-1)+1:npc*i,[1 3])=perfilNACA;
end
k=0;
for i=1:n 
    for j=1:epc
        k=k+1;
        T(epc*(i-1)+j,:)=ConectivitatsNACA(j,:)+npc*(i-1);
    end
end
nbc=length(conectivitatsEstructura); %Nombre de barres creuades
for i=2:n
    for j=1:nbc
        T((nbc)*(i-2)+j+k,:)=[conectivitatsEstructura(j,1)+npc*(i-2) conectivitatsEstructura(j,2)+npc*(i-1)];
    end
end

end

