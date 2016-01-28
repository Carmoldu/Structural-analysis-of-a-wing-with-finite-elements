function [ Fa0,Fa1,Fa2,Fa3 ] = matAero(xpoints,x0,v0,U,omega,env,c,n,alphaDeg )
    %% MATRIUS DE FORCES AERODINÀMIQUES
    %% Adaptació de la geometria
    
    for i=1:(n-1)*6
        for j=1:3
            X((i-1)*3+j)=xpoints(i+6,j);
        end
    end
    X=X';

    %% Definició de la matriu de forces aerodinàmiques


    b=c/2;        %m
    a=env/(n-1);
    
    rho=1.225;  %k/m^3
    %omega=10;   %frequència del mode principal
    k=omega*b/U; %freqüència reduida
    C=@(k) real(besselh(1,2,k)/(besselh(1,2,k)+1i*besselh(0,2,k))); 

    Lh=0;
    Lh1=2*pi*rho*U*b*C(k);
    Lh2=pi*rho*b^2;
    La=Lh1*U;
    La1=Lh2*U+Lh1*b*(1/2-a);
    La2=-Lh2*b*a;

    Mh=0;
    Mh1=2*U*(a+1/2)*C(k);
    Mh2=Lh2*b*a;
    Ma=Mh1*U;
    Ma1=-Lh2*U*b*(1/2-a)+Mh1*b*(1/2-a);
    Ma2=-Lh2*b^2*(1/8+a^2);

    %% Matriu de forces Aerodinàmiques
    LM=[Lh La Lh1 La1 Lh2 La2
        Mh Ma Mh1 Ma1 Mh2 Ma2];

    %% Transormació geomètrica


    zca=X(3)+(X(6*3)-X(3))*1/4; %Alltura de intesecció de la línia de sustentació nula amb c/4

    % Coeficients de h
    kh=zeros(1,6);
    kh(2)=(X(3*3)-zca)/(X(3*3)-X(3*2));
    kh(3)=-(X(2*3)-zca)/(X(3*3)-X(3*2));

    % Coeficients de alpha
    ka=zeros(1,6);
    ka(1)=1/3;
    ka(4)=1/6;
    ka(5)=1/6;
    ka(6)=1/3;

    % Distància a c/4
    cm(1)=-c/4;
    cm(2)=0;
    cm(3)=0;
    cm(4)=c/4;
    cm(5)=c/4;
    cm(6)=3/4*c;

    %% Generació de la matriu de transformació geomètrica
    Tx=zeros(2,6);
    Tx(1,:)=kh;
    for i=1:6
        if i~=2 && i~=3
            Tx(2,i)=ka(i)/cm(i);
        end
    end
    sumTx=sum(Tx(2,:));
    Tx(2,2)=-kh(2)*sumTx;
    Tx(2,3)=-kh(3)*sumTx;  

    %% Transformació de forces

    % Distribució de la sustentació
    r1=0.5;
    r2=3;
    funL =@(kL) [kL(1)+kL(2)+kL(3)+kL(4)+kL(5)+kL(6)-1
                kL(1)*cm(1)+kL(2)*cm(2)+kL(3)*cm(3)+kL(4)*cm(4)+kL(5)*cm(5)+kL(6)*cm(6)
                kL(2)-kL(3)
                kL(4)-kL(5)
                kL(1)-r1*kL(2)
                kL(5)-r2*kL(6)];
    kL=fsolve(funL,zeros(1,6));

    % Distribució del moment
    r3=3;
    funM =@(kM) [kM(1)+kM(2)+kM(3)+kM(4)+kM(5)+kM(6)-1
                kM(1)/cm(1)+kM(4)/cm(4)+kM(5)/cm(5)+kM(6)/cm(6)
                kM(2)
                kM(3)
                kM(4)-kM(5)
                kM(5)-r3*kM(6)];
    kM=fsolve(funM,zeros(1,6));

    %% Generació de la matriu de transformació de forces
    Tf=zeros(2,6);
    Tf(1,:)=-kL;
    for i=1:6
        if i~=2 && i~=3
            Tf(2,i)=kM(i)/cm(i);
        end
    end
    Tf=Tf';

    %% Concatenació de les matrius
    Txg=zeros(2,6*3);
    Tfg=zeros(6*3,2);

    for i=1:6
        Txg(:,i*3)=Tx(:,i);
        Tfg(i*3,:)=Tf(i,:);
    end

    TxG=zeros(4,(n-1)*6*3);
    TfG=zeros((n-1)*6*3,4);
    % for i=1:n-2
    %     TxG=[TxG Txg];
    %     TfG=[TfG;Tfg];
    % end
    dofc=6*3;
    var=2;
    for i=1:n-1
        TxG((i-1)*var+1:i*var,(i-1)*dofc+1:i*dofc)=Txg;
        TfG((i-1)*dofc+1:i*dofc,(i-1)*var+1:i*var)=Tfg;
        LM1((i-1)*var+1:i*var,(i-1)*var+1:i*var)=LM(:,1:2);
        LM2((i-1)*var+1:i*var,(i-1)*var+1:i*var)=LM(:,3:4);
        LM3((i-1)*var+1:i*var,(i-1)*var+1:i*var)=LM(:,5:6);
    end
    TX1=zeros((n-1)*4,(n-1)*6*3*2);
    TF1=zeros((n-1)*6*3*2,(n-1)*4);
    LMG1=zeros(var*(n-1)*2);

    TX1(1:end/2,1:end/2)=TxG;
    TX2=TX1;
    TX1(end/2+1:end,end/2+1:end)=TxG;

    LMG2=LMG1;
    LMG2(1:end/2,1:end/2)=LM3;

    LMG1(1:end/2,1:end/2)=LM1;
    LMG1(end/2+1:end,end/2+1:end)=LM2;


    TF1(1:end/2,1:end/2)=TfG;
    TF2=TF1;
    TF1(end/2+1:end,end/2+1:end)=TfG;


 
    alpha0=alphaDeg*pi/180; 
    alpha=zeros((n-1)*var,1);
    for i=1:n-1
        alpha(i*2)=alpha0;
    end

    Fa3=TfG*LM3*TxG;

    Fa2=TfG*LM2*TxG;

    Fa1=TfG*LM1*TxG;

    Fa0=TfG*LM1*(TxG*X+alpha)+TfG*LM2*TxG*v0;


end

