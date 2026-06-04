
clear;clc;

% INPUT FILE PATH
INP_FILE = 'SAMPLE_INPUT_FILE2D.yaml';

% READ THE INPUT FILE AND EXTRACT REQUIRED VARIABLES
INP_YAML = ReadYaml(INP_FILE);

X0 = INP_YAML.X0; % X COORDINATE OF THE BOTTOM-LEFT CORNER OF THE DOMAIN
Y0 = INP_YAML.Y0; % Y COORDINATE OF THE BOTTOM-LEFT CORNER OF THE DOMAIN
DX = cell2mat(INP_YAML.DX); % ELEMENT LENGTHS ALONG X
DY = cell2mat(INP_YAML.DY); % ELEMENT LENGTHS ALONG Y
NPE = INP_YAML.NPE; % NODES PER ELEMENT

%READING TIME PARAMETERS
ITEM = INP_YAML.ITEM;
DT = INP_YAML.TIMESTEP;
NTIME = INP_YAML.NTIME;

%READING STABILITY PARAMETERS
ALPHA =INP_YAML.ALPHA;
NEW = INP_YAML.NEW;

POINT = INP_YAML.POINT;

% PDE COEFFICIENTS...
A11 = INP_YAML.A11;
A22 = INP_YAML.A22;
A0 = INP_YAML.A0;
C1 = INP_YAML.C1;
C2 = INP_YAML.C2;
F  = INP_YAML.F;
[DATA.A0] = deal(A0{:});
[DATA.A10, DATA.A1X, DATA.A1Y, DATA.A1U, DATA.A1UX, DATA.A1UY] = ...
    deal(A11{:});
[DATA.A20, DATA.A2X, DATA.A2Y, DATA.A2U, DATA.A2UX, DATA.A2UY] = ...
    deal(A22{:});
[DATA.C10, DATA.C1X1, DATA.C1Y1] = deal(C1{:});
[DATA.C20, DATA.C2X1, DATA.C2Y1] = deal(C2{:});
[DATA.FX0, DATA.FX1, DATA.FY1] = deal(F{:});

DATA.A1 = ALPHA*DT;
DATA.A2 = (1-ALPHA)*DT;

if ITEM>1
    DATA.A3 = 2/(NEW*DT^2);
    DATA.A4 = DATA.A3*DT;
    DATA.A5 = 1/NEW - 1;
    DATA.A6 = (2*ALPHA)/(NEW*DT);
    DATA.A7 = (2*ALPHA)/NEW - 1;
    DATA.A8 = DT*((ALPHA/NEW)-1);
else
    DATA.A3 = 0;
    DATA.A4 = 0;
    DATA.A5 = 0;
    DATA.A6 = 0;
    DATA.A7 = 0;
    DATA.A8 = 0;
end

NGP  = INP_YAML.NGP;  % NUMBER OF GAUSS POINTS PER DIRECTION
NSPV = INP_YAML.NSPV; % NUMBER OF DIRICHLET BOUNDARY CONDITIONS
NSSV = INP_YAML.NSSV; % NUMBER OF NEUMANN BOUNDARY CONDITIONS
NSMB = INP_YAML.NSMB; % NUMBER OF ROBIN BOUNDARY CONDITIONS

% NODES AND DEGREES OF FREEDOM WHERE THE 
% PRIMARY VARIABLE IS PRESCRIBED...
ISPV(:,1) = cell2mat(INP_YAML.ISPV.NODES);
ISPV(:,2) = cell2mat(INP_YAML.ISPV.DOFS);
VSPV      = cell2mat(INP_YAML.VSPV);

% NODES AND DEGREES OF FREEDOM WHERE THE 
% SECONDARY VARIABLE IS PRESCRIBED...
if(NSSV>0)
    ISSV(:,1) = cell2mat(INP_YAML.ISSV.NODES);
    ISSV(:,2) = cell2mat(INP_YAML.ISSV.DOFS);
    VSSV      = cell2mat(INP_YAML.VSSV);
end

% NODES AND DOFS WITH PRESCRIBED 
% ROBIN (MIXED) BOUNDARY CONDITIONS...
if(NSMB>0)
    ISMB(:,1) = cell2mat(INP_YAML.ISMB.NODES);
    ISMB(:,2) = cell2mat(INP_YAML.ISMB.DOFS);
    UREF      = cell2mat(INP_YAML.UREF);
    BETA0     = cell2mat(INP_YAML.BETA0);
    BETAU     = cell2mat(INP_YAML.BETAU);
end

% NONLINEAR ANALYSIS PARAMETERS
NONLIN  = INP_YAML.NONLIN;  % SOLVER MODE (0-LINEAR; 1-PICARD; 2-NEWTON)
ITERMAX = INP_YAML.ITERMAX; % MAXIMUM NO. OF ITERATIONS BEFORE TERMINATION
EPSILON = INP_YAML.EPSILON; % CONVERGENCE TOLERANCE
GAMMA   = INP_YAML.GAMMA;   % RELAXATION PARAMETER


% M E S H   G E N E R A T I O N
% .........................................................................

% CALL THE 2D MESHING FUNCTION
NX = numel(DX);
NY = numel(DY);
[NOD, GLXY, NEM, NNM] = MESH2DR(NX,NY,NPE,DX,DY,X0,Y0);

% REMARKS:
% GLXY - ARRAY OF NODE COORDINATES
% NOD  - ELEMENT CONNECTIVITY MATRIX
% NNM  - TOTAL NUMBER OF NODES IN THE MESH
% NPE  - NUMBER OF NODES PER ELEMENT

% PRECOMPUTATION OF SHAPE FUNCTIONS AND THEIR 
% DERIVATIVES AT INTEGRATION POINTS
GFILE = "GAUSS.g"; % DATABASE OF GAUSS QUADRATURE DATA
[GAUSWT,SFL_ARRAY,DSFL_ARRAY] = PRECOMPUTE_SF2D(GFILE,NGP,NPE);


% F I N I T E   E L E M E N T   S O L V E R
% .........................................................................

NDF = 1;       % NUMBER OF DEGREES OF FREEDOM PER NODE
NEQ = NNM*NDF; % NUMBER OF EQUATIONS IN THE GLOBAL SYSTEM
NET = NPE*NDF; % NUMBER OF ELEMENT LEVEL EQUATIONS

% INITIALIZE THE SOLUTION ARRAYS
GCU = zeros(NEQ,1);   % GLOBAL SOLUTION ARRAY (CURRENT ITERATION)
GPU = zeros(NEQ,1);   % GLOBAL SOLUTION ARRAY (PREVIOUS ITERATION)

GLPU = zeros(NEQ,1); %PREVIOUS TIMESTEP SOLUTION
GLPV = zeros(NEQ,1);
GLPA = zeros(NEQ,1);
GLCU = zeros(NEQ,1); %CURRENT TIMESTEP SOLUTION
GLCV = zeros(NEQ,1);
GLCA = zeros(NEQ,1);
GLC = zeros(NEQ,1); %CURRENT TIMESTEP & CURRENT ITERATION SOLUTION

TIME_HIST = zeros(NTIME,1);

%TIME LOOP
for T = 1:NTIME
    TIME = T*DT;
    
    % NONLINEAR SOLUTION ITERATION LOOP
    for ITER = 1:ITERMAX
        
        % INITIALIZE GLOBAL STIFFNESS AND FORCE...
        GLK = zeros(NEQ,NEQ); % GLOBAL 'STIFFNESS' MATRIX
        GLF = zeros(NEQ,1);   % GLOBAL 'FORCE' VECTOR
    
        % LOOP OVER ELEMENTS
        for N=1:NEM
            NODES = NOD(N,:);      % NODES OF CURRENT ELEMENT
            ELXY  = GLXY(NODES,:); % NODAL COORDINATES OF CURRENT ELEMENT
        
            % CONSTRUCT GLOBAL-TO-LOCAL DOF MAP
            INDXS = zeros(1,NET);
            S = 1;
            for ND = 1:NDF
                INDXS(S:NDF:end) = NDF*NODES - (NDF - S)*ones(1,NPE);
                S = S + 1;
            end
            
            % Compute the element-level solution arrays
            ELPU = GLPU(INDXS);
            ELPV = GLPV(INDXS);
            ELPA = GLPA(INDXS);
            ELCU = GLCU(INDXS);
            ELCV = GLCV(INDXS);
            ELCA = GLCA(INDXS);

            % ELEMENT SOLUTION ARRAY AFTER APPLYING RELAXATION
            %ELU = (1-GAMMA)*GPU(INDXS) + GAMMA*GCU(INDXS);
            
            % CALL ELEMENT-LEVEL ROUTINE TO COMPUTE 
            % THE ELEMENT COEFFICIENT ARRAYS
            [ELK,ELF,ELK_EFF,ELF_EFF] = ELEMAT2D(DATA,SFL_ARRAY,DSFL_ARRAY,GAUSWT,NONLIN,...
                                  ELXY,ELPU,ELPV,ELPA,ELCU,ELCV,ELCA,ITEM);
    
            % ASSEMBLE ELEMENT-LEVEL EQUATIONS INTO THE GLOBAL SYSTEM
            GLK(INDXS,INDXS) = GLK(INDXS,INDXS)+ELK_EFF;
            GLF(INDXS) = GLF(INDXS)+ELF_EFF;
        end
    
        % MODIFY GLOBAL MATRICES (GLK, GLF) TO 
        % ENFORCE DIRICHLET BOUNDARY CONDITIONS
        if(NSPV>0)
    
            for NP=1:NSPV
                NB = (ISPV(NP,1)-1)*NDF+ISPV(NP,2);
                for J=1:NEQ
                    GLK(NB,J)  = 0;
                    GLK(NB,NB) = 1;
                    GLF(NB)    = VSPV(NP);
                end
            end
        end
        
        % MODIFY THE FORCE VECTOR TO INCLUDE SPECIFIED NONZERO 
        % SECONDARY VARIABLES IN THE FULL MATRIX SYSTEM
        if(NSSV>0)
            for NS=1:NSSV
                NB = (ISSV(NS,1)-1)*NDF+ISSV(NS,2);
                GLF(NB) = GLF(NB) + VSSV(NS);
            end
        end
    
        % IMPLEMENT THE SPECIFIED ROBIN BOUNDARY CONDITIONS
        if(NSMB>0)
            for MB=1:NSMB
                NB = (ISMB(MB,1)-1)*NDF+ISMB(MB,2);
                if(NONLIN<=1)
                    GLK(NB,NB)=GLK(NB,NB)+BETA0(MB)+BETAU(MB)*GCU(NB);
                    GLF(NB)=GLF(NB)+UREF(MB)*(BETA0(MB)+BETAU(MB)*GCU(NB));
                else
                    GLK(NB,NB)= GLK(NB,NB)+BETA0(MB)+2.0*BETAU(MB)*GCU(NB)...
     	                                      -UREF(MB)*BETAU(MB);
                    GLF(NB)= GLF(NB)+UREF(MB)*(BETA0(MB) ...
                                +BETAU(MB)*GCU(NB))-(BETA0(MB) ...
                                +BETAU(MB)*GCU(NB))*GCU(NB);
                end
            end
        end
    
        % SOLVE SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS. 
        SOLU = GLK\GLF;
        %GLPU=GLCU;
        if(NONLIN==0)
            GLCU=SOLU;
            disp("LINEAR SOLUTION COMPLETE.");
            break;
        end
        if NONLIN==1
            GPU = GLCU;
            GLCU = SOLU;
            DU = GLCU - GPU;
            ERR = norm(DU)/norm(GLCU);
            if ERR<EPSILON
                break;
            end
        elseif NONLIN>1
            GPU = GLCU;
            GLCU = GLCU + SOLU;
            if norm(SOLU) < EPSILON
                break;
            end
        end
       
    end
    %GPU = GLPU;
    % GLPU = GLCU;

    if ITEM>1
        GLCA = DATA.A3 *(GLCU-GLPU) - DATA.A4*GLPV - DATA.A5*GLPA;
        GLCV = GLPV +DATA.A1*GLCA + DATA.A2*GLPA;
        GLPA = GLCA;
        GLPV = GLCV;
    end
    GLPU = GLCU;
    disp(GLCU);
    TIME_HIST(T,1) = GLCU(POINT);
end
% END OF THE NONLINEAR SOLUTION ITERATION LOOP

function [ELK,ELF,ELK_EFF,ELF_EFF] = ELEMAT2D(DATA,SFL_ARRAY,DSFL_ARRAY,GAUSWT,NONLIN,...
                              ELXY,ELPU,ELPV,ELPA,ELCU,ELCV,ELCA,ITEM)
    
    NDF = 1;            % DEGREES OF FREEDOM PER NODE
    NPE = size(ELXY,1); % NODES PER ELEMENT
    NET = NPE*NDF;      % NUMBER OF ELEMENT LEVEL EQUATIONS
    
    % EXTRACT PDE COEFFICIENTS FROM THE DATA STRUCT
    [A0] = deal(DATA.A0);
    [A10, A1X, A1Y, A1U, A1UX, A1UY] = ... 
        deal(DATA.A10, DATA.A1X, DATA.A1Y, DATA.A1U, DATA.A1UX, DATA.A1UY);
    [A20, A2X, A2Y, A2U, A2UX, A2UY] = ... 
        deal(DATA.A20, DATA.A2X, DATA.A2Y, DATA.A2U, DATA.A2UX, DATA.A2UY);
    [FX0, FX1, FY1] = ...
        deal(DATA.FX0, DATA.FX1, DATA.FY1);
    [C10, C1X1, C1Y1] = deal(DATA.C10, DATA.C1X1, DATA.C1Y1);
    [C20, C2X1, C2Y1] = deal(DATA.C20, DATA.C2X1, DATA.C2Y1);
    
    % INITIALIZE ELEMENT COEFFICIENTS
    ELK = zeros(NET,NET);
    ELF = zeros(NET,1);
    ELM = zeros(NET, NET);
    ELC = zeros(NET, NET);

    if(NONLIN>1)
        % PART OF THE ELEMENT-LEVEL TANGENT MATRIX THAT ADDS TO ELK
        TAN = zeros(NPE,NPE);
    end
    
    % NUMBER OF GAUSS POINTS
    NGP2 = size(GAUSWT,1);

    % INTEGRATION LOOP
    for NG = 1:NGP2
        
        % SHAPE FUNCTION DATA AT CURRENT GAUSS POINT
        SFL   = SFL_ARRAY(:,NG);     % SHAPE FUNCTIONS
        DSFL  = DSFL_ARRAY(:,:,NG);  % SHAPE FUNCTION DERIVATIVES (LOCAL)
        JACOB = ELXY'*DSFL;          % JACOBIAN OF TRANSFORMATION
        GDSFL = DSFL/JACOB;          % SHAPE FUNCTION DERIVATIVES (GLOBAL)
        CNST  = det(JACOB)*GAUSWT(NG);
        
        % GLOBAL COORDINATES OF THE CURRENT GAUSS POINT
        X = dot(ELXY(:,1),SFL);
        Y = dot(ELXY(:,2),SFL);

        S00 = SFL * SFL' * CNST;
        
        if(NONLIN>0)
            U   = dot(ELCU,SFL);        % 'u'
            DUX = dot(ELCU,GDSFL(:,1)); % 'du/dx'
            DUY = dot(ELCU,GDSFL(:,2)); % 'du/dy'
        end
        
        A11 = A10 + A1X*X + A1Y*Y;
        A22 = A20 + A2X*X + A2Y*Y;
        if ITEM>1
            C1 = C10 + C1X1*X + C1Y1*Y;
            C2 = C20 + C2X1*X + C2Y1*Y;
        else
            C1 = C10 + C1X1*X + C1Y1*Y;
            C2 = 0;
        end
        % C1 = C10 + C1X1*X + C1Y1*Y;
        % C2 = C20 + C2X1*X + C2Y1*Y;
        
        F = FX0 + FX1*X + FY1*Y;
        if(NONLIN>0)
           A11 = A10 + A1X*X + A1Y*Y + A1U*U + A1UX*DUX + A1UY*DUY;
           A22 = A20 + A2X*X + A2Y*Y + A2U*U + A2UX*DUX + A2UY*DUY; 
        end
        
        % COMPUTE ELEMENT STIFFNESS, MASS, DAMPING MATRIX AND FORCE VECTOR

        if ITEM>1
            ELM = ELM + C2*S00;
        else
            ELM = 0;
        end

        ELC = ELC + C1*S00;
        ELK = ELK + A0*S00 +...
                    A11*(GDSFL(:,1))*GDSFL(:,1)'*CNST + ...
                    A22*(GDSFL(:,2))*GDSFL(:,2)'*CNST;
        ELF = ELF + F*SFL*CNST;
        
        if(NONLIN>1)

            TAN = TAN + DUX*CNST*(A1U*GDSFL(:,1)*SFL' + A1UX*(GDSFL(:,1))*GDSFL(:,1)' + A1UY*(GDSFL(:,1))*GDSFL(:,2)')...
                + DUY*CNST*(A2U*GDSFL(:,2)*SFL' + A1UX*(GDSFL(:,2))*GDSFL(:,1)' + A1UY*(GDSFL(:,2))*GDSFL(:,2)');

        end
    end
    if ITEM>1
        AS = DATA.A3*ELPU + DATA.A4*ELPV + DATA.A5*ELPA;   % A_n
        BS = DATA.A6*ELPU + DATA.A7*ELPV + DATA.A8*ELPA;   % B_n
        ELK_EFF = ELK + DATA.A3*ELM + DATA.A6*ELC;
        ELF_EFF = ELF + ELM*AS +ELC*BS;
        if NONLIN > 1
            ELK_EFF = ELK_EFF + TAN;           % tangent
            ELF_EFF = ELF_EFF - ELK_EFF*ELCU;     % = -residual
        end
    else
        ELK_EFF = ELC + DATA.A1*ELK;
        ELF_EFF = (ELC - DATA.A2*ELK)*ELPU + (DATA.A1 + DATA.A2)*ELF;
    end
    % ELK_EFF = ELK + DATA.A3*ELM + DATA.A6*ELC;
    % ELF_EFF = ELF + ELM*AS +ELC*BS;
    
    % FOR NEWTON'S METHOD
    % if(NONLIN>1)
    %     % ELEMENT-LEVEL RESIDUAL
    %     ELF = ELF - ELK*ELU;
    %     % ELEMENT-LEVEL TANGENT MATRIX
    %     ELK = ELK + TAN;
    % end
end


function [NOD, GLXY, NEM, NNM] = MESH2DR(NX,NY,NPE,DX,DY,X0,Y0)
    if(NPE==4)
        IEL = 1;
    else
        IEL = 2;
    end
    NEX1 = NX+1;
    NEY1 = NY+1;
    NXX  = IEL*NX;
    NYY  = IEL*NY;
    NXX1 = NXX + 1;
    NYY1 = NYY + 1;
    NEM  = NX*NY;
    NNM  = NXX1*NYY1;
    if(NPE==8)
        NNM = NXX1*NYY1 - NX*NY;
    end
    GLXY = zeros(NNM,2);
    K0 = 0;
    if (NPE==9) 
        K0=1;
    end
    NOD(1,1) = 1;
    NOD(1,2) = IEL+1;
    NOD(1,3) = NXX1+(IEL-1)*NEX1+IEL+1;
    if(NPE==9) 
        NOD(1,3)=4*NX+5;
    end
    NOD(1,4) = NOD(1,3) - IEL;
    if(NPE>4)
        NOD(1,5) = 2;
        NOD(1,6) = NXX1 + (NPE-6);
        NOD(1,7) = NOD(1,3) - 1;
        NOD(1,8) = NXX1+1;
        if(NPE==9)
            NOD(1,9)=NXX1+2;
        end
    end
    if(NY>1) 
        M = 1;
        for N = 2:NY
            L = (N-1)*NX + 1;
            for I = 1:NPE
                NOD(L,I) = NOD(M,I)+NXX1+(IEL-1)*NEX1+K0*NX;
            end
            M=L;
        end
    end
    if(NX>1)
        for NI = 2:NX
            for I = 1:NPE
                K1 = IEL;
                if(I==6 || I==8)
                    K1=1+K0;
                end
                NOD(NI,I) = NOD(NI-1,I)+K1;
            end
            M = NI;
            for NJ = 2:NY
                L = (NJ-1)*NX+NI;
                for J = 1:NPE
                    NOD(L,J) = NOD(M,J)+NXX1+(IEL-1)*NEX1+K0*NX;
                end
                M = L;
            end
        end
    end
  
    DX(NEX1) = 0.0;
    DY(NEY1) = 0.0;
    XC = X0;
    YC=Y0;
    if(NPE==8)
        for NI = 1:NEY1 
            I = (NXX1+NEX1)*(NI-1)+1;
            J = 2*NI-1;
            GLXY(I,1) = XC;
            GLXY(I,2) = YC;
            for NJ = 1:NX
                DELX=0.5*DX(NJ);
                I=I+1;
                GLXY(I,1) = GLXY(I-1,1)+DELX;
                GLXY(I,2) = YC;
                I=I+1;
                GLXY(I,1) = GLXY(I-1,1)+DELX;
                GLXY(I,2) = YC;
            end
            if(NI<=NY)
                I = I+1;
                YC= YC+0.5*DY(NI);
                GLXY(I,1) = XC;
                GLXY(I,2) = YC;
                for  II = 1:NX
                    I = I+1;
                    GLXY(I,1) = GLXY(I-1,1)+DX(II);
                    GLXY(I,2) = YC;
                end
            end
            YC = YC+0.5*DY(NI); 
        end
    else
        YC=Y0;
        for NI = 1:NEY1
            XC = X0;
            I = NXX1*IEL*(NI-1);
            for NJ = 1:NEX1
                I=I+1;
                GLXY(I,1) = XC;
                GLXY(I,2) = YC;
                if(NJ<NEX1)
                    if(IEL==2)
                        I=I+1;
                        XC = XC + 0.5*DX(NJ);
                        GLXY(I,1) = XC;
                        GLXY(I,2) = YC;
                    end
                end
                XC = XC + DX(NJ)/IEL;
            end
            XC = X0;
            if(IEL==2)
                YC = YC + 0.5*DY(NI);
                for NJ = 1:NEX1
                    I=I+1;
                    GLXY(I,1) = XC;
                    GLXY(I,2) = YC;
                    if(NJ<NEX1)
                        I=I+1;
                        XC = XC + 0.5*DX(NJ);
                        GLXY(I,1) = XC;
                        GLXY(I,2) = YC;
                    end
                    XC = XC + 0.5*DX(NJ);
                end
            end
            YC = YC + DY(NI)/IEL;
        end
    end
end

function [SFL, DSFL] = INTERPLN2D(NPE,XI,ETA)
    
    NP = [1,2,3,4,5,7,6,8,9];
    XNODE = [-1.0D0, 1.0D0,1.0D0, -1.0D0, 0.0D0, 1.0D0, 0.0D0, -1.0D0,...
             0.0D0; -1.0D0,-1.0D0, 1.0D0,1.0D0, -1.0D0, 0.0D0, 1.0D0,...
             0.0D0,0.0D0]';
    
    SFL = zeros(NPE,1);
    DSFL = zeros(2,NPE);
    if(NPE==4)
    % LINEAR LAGRANGE INTERPOLATION FUNCTIONS FOR FOUR-NODE ELEMENT
        for I = 1:NPE
            XP  = XNODE(I,1);
            YP  = XNODE(I,2);
            XI0 = 1.0+XI*XP;
            ETA0=1.0+ETA*YP;
            SFL(I)   = 0.25*FNC(XI0,ETA0);
            DSFL(1,I)= 0.25*FNC(XP,ETA0);
            DSFL(2,I)= 0.25*FNC(YP,XI0);
         end
    elseif(NPE==8)
    % QUADRATIC LAGRANGE INTERPOLATION FUNCTIONS FOR EIGHT-NODE ELEMENT
        for I = 1:NPE
            NI   = NP(I);
            XP   = XNODE(NI,1);
            YP   = XNODE(NI,2);
            XI0  = 1.0+XI*XP;
            ETA0 = 1.0+ETA*YP;
            XI1  = 1.0-XI*XI;
            ETA1 = 1.0-ETA*ETA;
            if(I<=4)
                SFL(NI)    = 0.25*FNC(XI0,ETA0)*(XI*XP+ETA*YP-1.0);
                DSFL(1,NI) = 0.25*FNC(ETA0,XP)*(2.0*XI*XP+ETA*YP);
                DSFL(2,NI) = 0.25*FNC(XI0,YP)*(2.0*ETA*YP+XI*XP);
            else
                if(I<=6)
                    SFL(NI)    = 0.5*FNC(XI1,ETA0);
                    DSFL(1,NI) = -FNC(XI,ETA0);
                    DSFL(2,NI) = 0.5*FNC(YP,XI1);
                else
                    SFL(NI)    = 0.5*FNC(ETA1,XI0);
                    DSFL(1,NI) = 0.5*FNC(XP,ETA1);
                    DSFL(2,NI) = -FNC(ETA,XI0);
                end
            end
        end
    elseif(NPE==9)
    % QUADRATIC LAGRANGE INTERPOLATION FUNCTIONS FOR NINE-NODE ELEMENT
        for I = 1:NPE
            NI   = NP(I);
            XP   = XNODE(NI,1);
            YP   = XNODE(NI,2);
            XI0  = 1.0+XI*XP;
            ETA0 = 1.0+ETA*YP;
            XI1  = 1.0-XI*XI;
            ETA1 = 1.0-ETA*ETA;
            XI2  = XP*XI;
            ETA2 = YP*ETA;
            if(I <=4)
               SFL(NI)   = 0.25*FNC(XI0,ETA0)*XI2*ETA2;
               DSFL(1,NI)= 0.25*XP*FNC(ETA2,ETA0)*(1.0+2.0*XI2);
               DSFL(2,NI)= 0.25*YP*FNC(XI2,XI0)*(1.0+2.0*ETA2);
            else
               if(I<=6)
                  SFL(NI)    = 0.5*FNC(XI1,ETA0)*ETA2;
                  DSFL(1,NI) = -XI*FNC(ETA2,ETA0);
                  DSFL(2,NI) = 0.5*FNC(XI1,YP)*(1.0+2.0*ETA2);
               else
                  if(I<=8)
                     SFL(NI)    = 0.5*FNC(ETA1,XI0)*XI2;
                     DSFL(2,NI) = -ETA*FNC(XI2,XI0);
                     DSFL(1,NI) = 0.5*FNC(ETA1,XP)*(1.0+2.0*XI2);
                  else
                     SFL(NI)    = FNC(XI1,ETA1);
                     DSFL(1,NI) = -2.0*XI*ETA1;
                     DSFL(2,NI) = -2.0*ETA*XI1;
                  end
               end
            end
        end
    end
end

function Y = FNC(A,B)
    Y = A*B;
end



function [GAUSWT,SFL_ARRAY,DSFL_ARRAY] = PRECOMPUTE_SF2D(GFILE,NGP,NPE)
    [GAUSPT,GAUSWT] = GAUSS(GFILE,NGP);
    NGP2 = NGP*NGP;
    GAUSWT = reshape(GAUSWT*GAUSWT',[NGP2,1]);
    [IPET,IPXI] = meshgrid(GAUSPT,GAUSPT);
    IPTS = [IPXI(:),IPET(:)];
    SFL_ARRAY = zeros(NPE,NGP2);
    DSFL_ARRAY = zeros(NPE,2,NGP2);
    for I = 1:NGP2
        XI = IPTS(I,1);
        ETA = IPTS(I,2);
        [SFL, DSFL] = INTERPLN2D(NPE,XI,ETA);
        SFL_ARRAY(:,I) = SFL;
        DSFL_ARRAY(:,:,I) = transpose(DSFL);
    end
end

function [GAUSPT,GAUSWT] = GAUSS(GFILE,NGP)
    FILEID  = fopen(GFILE);
    TLINE   = fgetl(FILEID);
    TLINES  = cell(50,1);
    K = 1;
    while ischar(TLINE)
        TLINES{K,1} = TLINE; TLINE = fgetl(FILEID); K = K + 1;
    end
    fclose(FILEID);
    LINE = TLINES{2};
    CONT = textscan(LINE,'%f','Delimiter',',');
    SF_S = CONT{1};
    SF   = SF_S(NGP);
    GAUSPT = zeros(NGP,1);
    GAUSWT = zeros(NGP,1);
    for NI = 1:NGP
        LINE = TLINES{SF};
        CONT = textscan(LINE,'%f %f %f');
        GAUSWT(NI) = CONT{2};
        GAUSPT(NI) = CONT{3};
        SF = SF + 1;
    end
end

