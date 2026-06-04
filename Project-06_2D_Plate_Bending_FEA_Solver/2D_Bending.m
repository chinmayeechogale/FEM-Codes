
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
Q    = INP_YAML.Q;


NGPF  = INP_YAML.NGP; % NUMBER OF GAUSS POINTS
NGPR  = INP_YAML.LGP; %Number of Gauss points for reduced integration

E1  = INP_YAML.E1;     % Young's Modulus
E2  = INP_YAML.E2;  
G12 = INP_YAML.G12;     
G13 = INP_YAML.G13;
G23 = INP_YAML.G23;
H = INP_YAML.H;
Ks   = INP_YAML.Ks;    % Shear Correction Coefficient
Nu12 = INP_YAML.Nu12;    % Poisson's Ratio

Nu21 = Nu12*E2/E1;

%Calculating coefficients
DATA.H = H;
DATA.A11 = (E1*H)/(1-Nu12*Nu21);
DATA.A12 = Nu21*DATA.A11;
DATA.A22 = DATA.A11*(E2/E1);
DATA.D11 = DATA.A11*(H*H/12);
DATA.D12 = Nu21*DATA.D11;
DATA.D22 = DATA.D11*(E2/E1);
DATA.A66 = G12*H;
DATA.A44 = Ks*G23*H;
DATA.A55 = Ks*G13*H;
DATA.D66 = G12*(H^3/12);

DATA.Q11 = E1/(1-(Nu12*Nu21));
DATA.Q12 = Nu12 *E2/(1-(Nu12*Nu21));
DATA.Q21 = DATA.Q12;
DATA.Q22 = E2/(1-(Nu12*Nu21));
DATA.Q44 = G12;
DATA.Q55 = G12;
DATA.Q66 = G12;

% DATA.A11 = 1.0667e4;
% DATA.A12 = 0.2667e4;
% DATA.A22 = 1.0667e4;
% DATA.D11 = 8.8889e-2;
% DATA.D12 = 2.2222e-2;
% DATA.D22 = 8.8889e-2;
% DATA.A44 = 0.4e4*Ks;
% DATA.A55 = 0.4e4*Ks;
% DATA.A66 = 0.4e4;
% DATA.D66 = 3.3333e-2;
% 
% % PDE COEFFICIENTS...
% A1 = INP_YAML.A1;
% A2 = INP_YAML.A2;
FX  = INP_YAML.FX;
FY = INP_YAML.FY;

[DATA.FX0, DATA.FX1, DATA.FY1] = deal(FX{:});
[DATA.FY0, DATA.FX1, DATA.FY1] = deal(FY{:});
[QZ0, QZX1, QZY1] = deal(Q{:});
DATA.QZ = @(X,Y) QZ0+QZX1*X+QZY1*Y;
DATA.FX = @(X,Y) FX0 + FX1*X +FY1*Y; % << FIX THIS LATER IF THERE IS AXIAL DISTRIBUTED LOAD
DATA.FY = @(X,Y) FY0 + FY1*X +FY1*Y;

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
LSTEP   = INP_YAML.LSTEP;
DP = cell2mat(INP_YAML.DP);
RECDOF = INP_YAML.RECDOF;


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
[FULL.GAUSWT,FULL.SFL,FULL.DSFL] = PRECOMPUTE_SF2D(GFILE,NGPF,NPE);
[REDU.GAUSWT,REDU.SFL,REDU.DSFL] = PRECOMPUTE_SF2D(GFILE,NGPR,NPE);


%[FULL.SFL,FULL.DSFL,FULL.GAUSWT] = PRECOMPUTE_SF(GFILE,NGPF,NPE);

%[REDU.SFL,REDU.DSFL,REDU.GAUSWT] = PRECOMPUTE_SF(GFILE,NGPR,NPE);


% F I N I T E   E L E M E N T   S O L V E R
% .........................................................................

NDF = 5;       % NUMBER OF DEGREES OF FREEDOM PER NODE
NEQ = NNM*NDF; % NUMBER OF EQUATIONS IN THE GLOBAL SYSTEM
NET = NPE*NDF; % NUMBER OF ELEMENT LEVEL EQUATIONS

% INITIALIZE THE SOLUTION ARRAYS
GCU = zeros(NEQ,1);   % GLOBAL SOLUTION ARRAY (CURRENT ITERATION)
GPU = zeros(NEQ,1);   % GLOBAL SOLUTION ARRAY (PREVIOUS ITERATION)
W_MAX_ARRAY = [];
SXX_MAXX_ARRAY = [];
SXY_MAX_ARRAY = [];

ITER = 0; % ITERATION COUNTER
P = 0;
Q0 = E2*H^4/10^4;
for NL = 1:LSTEP
    CONVG = false;
    if(NONLIN)
        P = P + DP(NL);
    else
        P = 1;
    end
    F0 = P*Q0;
% NONLINEAR SOLUTION ITERATION LOOP
    for ITER = 1:ITERMAX
        
        % INITIALIZE GLOBAL STIFFNESS AND FORCE...
        GLK = zeros(NEQ,NEQ); % GLOBAL 'STIFFNESS' MATRIX
        GLF = zeros(NEQ,1);   % GLOBAL 'FORCE' VECTOR
        
        % STRESS_TOP_ALL = zeros(NEM,size(REDU.GAUSWT,1),5);
        % STRESS_BOT_ALL = zeros(NEM,size(REDU.GAUSWT,1),5);
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
    
            % ELEMENT SOLUTION ARRAY AFTER APPLYING RELAXATION
            ELU = (1-GAMMA)*GPU(INDXS) + GAMMA*GCU(INDXS);
            

            % CALL ELEMENT-LEVEL ROUTINE TO COMPUTE 
            % THE ELEMENT COEFFICIENT ARRAYS
            [ELK,ELF, STRS_TOP, STRS_BOT] = ELEMAT2D_FDST(DATA,ELU,ELXY,NPE, ...
                                  FULL,REDU,NGPF,NGPR, ...
                                  F0,NONLIN);
            % STRESS_TOP_ALL(N,:,:) = STRS_TOP;
            % STRESS_BOT_ALL(N,:,:) = STRS_BOT;
            % 
            % SXX_MAX = max(abs([STRESS_TOP_ALL(:,1);STRESS_BOT_ALL(:,1)]));
            % SXY_MAX = max(abs([STRESS_TOP_ALL(:,3);STRESS_BOT_ALL(:,3)]));
            % 
            % SXX_MAXX_ARRAY(1+end) = SXX_MAX;
            % SXY_MAX_ARRAY(1+end) = SXY_MAX;

            % ASSEMBLE ELEMENT-LEVEL EQUATIONS INTO THE GLOBAL SYSTEM
            GLK(INDXS,INDXS) = GLK(INDXS,INDXS)+ELK;
            GLF(INDXS) = GLF(INDXS)+ELF;
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
    
        GPU=GCU;
        if(NONLIN==0)
            GCU=SOLU;
            disp("LINEAR SOLUTION COMPLETE.");
            break;
        end
        if NONLIN==1
            GPU = GCU;
            GCU = SOLU;
            % DU = GCU - GPU;
            % ERR = DU/GCU;
            % if ERR<EPSILON
            %     break;
            % end
        elseif NONLIN>1
            GPU = GCU;
            GCU = GCU + SOLU;
            % if norm(SOLU) < EPSILON
            %     break;
            % end
        end
        DELTAU = GCU(:,1) - GPU(:,1);
        NUMER  = dot(DELTAU,DELTAU);
        DENOM  = dot(GCU,GCU);
        ERROR  = sqrt(NUMER/DENOM);
        if(ERROR<=EPSILON)
            CONVG = true;
            disp(NL);
            disp("LOAD STEP CONVERGED");
            disp(ITER);
            W = GCU(RECDOF);
            W_MAX_ARRAY(end+1) = W;
            disp(GCU(RECDOF));

            %disp(SGXX_BAR);
            %disp(W_MAX_ARRAY);
            break;
        end   
    end
    STRESS_TOP_ALL = zeros(NEM,size(REDU.GAUSWT,1),5);
    STRESS_BOT_ALL = zeros(NEM,size(REDU.GAUSWT,1),5);

    for N = 1:NEM
        NODES = NOD(N,:);      % NODES OF CURRENT ELEMENT
        ELXY  = GLXY(NODES,:); % NODAL COORDINATES OF CURRENT ELEMENT
        
            % CONSTRUCT GLOBAL-TO-LOCAL DOF MAP
         INDXS = zeros(1,NET);
         S = 1;
         for ND = 1:NDF
            INDXS(S:NDF:end) = NDF*NODES - (NDF - S)*ones(1,NPE);
            S = S + 1;
         end
        ELU = GCU(INDXS);
        [~,~,STRS_TOP,STRS_BOT] = ELEMAT2D_FDST(DATA,ELU,ELXY,NPE,FULL,REDU,NGPF,NGPR,F0,NONLIN);
        STRESS_TOP_ALL(N,:,:) = STRS_TOP;
        STRESS_BOT_ALL(N,:,:) = STRS_BOT;
    end

    SXX_MAX = max(abs([STRESS_TOP_ALL(:,:,1); STRESS_BOT_ALL(:,:,1)]),[],'all');
    SXY_MAX = max(abs([STRESS_TOP_ALL(:,:,3); STRESS_BOT_ALL(:,:,3)]),[],'all');
    SXX_MAXX_ARRAY(1+end) = SXX_MAX;
    SXY_MAX_ARRAY(1+end) = SXY_MAX;

end
    
% END OF THE NONLINEAR SOLUTION ITERATION LOOP


function [ELK,ELF, STRS_TOP, STRS_BOT] = ELEMAT2D_FDST(DATA,ELU,ELXY,NPE, ...
                                  FULL,REDU,NGPF,NGPR, ...
                                  F0,NONLIN)
    A11 = DATA.A11;
    A12 = DATA.A12;
    A22 = DATA.A22;
    D11 = DATA.D11;
    D12 = DATA.D12;
    D22 = DATA.D22;
    A66 = DATA.A66;
    A44 = DATA.A44;
    A55 = DATA.A55;
    D66 = DATA.D66;
    H = DATA.H;

    [FX0,FX1,FY1]  = deal(DATA.FX0, DATA.FX1, DATA.FY1);
    %[QZ0,QZX1,QZY1]  = deal(DATA.QZ0, DATA.QZX1, DATA.QZY1);
    
    NDF = 5;
    NET = NPE*NDF;
    K = 0;

    ELK   = zeros(NET,NET);
    ELF   = zeros(NET,1);
    ELK11 = zeros(NPE,NPE);
    ELK12 = zeros(NPE,NPE);
    %ELK21 = zeros(NPE,NPE);
    ELK22 = zeros(NPE,NPE);
    ELK23 = zeros(NPE,NPE);
    ELK32 = zeros(NPE,NPE);
    ELK33 = zeros(NPE,NPE);
    ELK34 = zeros(NPE,NPE);
    ELK35 = zeros(NPE,NPE);
    ELK24 = zeros(NPE,NPE);
    ELK25 = zeros(NPE,NPE);
    ELK13 = zeros(NPE,NPE);
    ELK31 = zeros(NPE,NPE);
    ELK14 = zeros(NPE,NPE);
    ELK15 = zeros(NPE,NPE);
    ELK41 = zeros(NPE,NPE);
    ELK42 = zeros(NPE,NPE);
    %ELK43 = zeros(NPE,NPE);
    ELK44 = zeros(NPE,NPE);
    ELK45 = zeros(NPE,NPE);
    ELK51 = zeros(NPE,NPE);
    ELK52 = zeros(NPE,NPE);
    %ELK53 = zeros(NPE,NPE);
    %ELK54 = zeros(NPE,NPE);
    ELK55 = zeros(NPE,NPE);

    ELF1 = zeros(NPE,1);
    ELF2 = zeros(NPE,1);
    ELF3 = zeros(NPE,1);

    NST = size(REDU.GAUSWT,1);
    STRS_TOP = zeros(NST,5);   % [sxx syy txy txz tyz]
    STRS_BOT = zeros(NST,5);   % [sxx syy txy txz tyz]
    

    if(NONLIN>1)
        TAN   = zeros(NET,NET);
        TAN13 = zeros(NPE,NPE);
        TAN23 = zeros(NPE,NPE);
        TAN33 = zeros(NPE,NPE);
    end

    for NG = 1:size(FULL.GAUSWT,1)
        SFL    = FULL.SFL(:,NG);
        DSFL   = FULL.DSFL(:,:,NG);
        J = ELXY'*DSFL;
        DETJ = det(J);
        GDSFL  = DSFL/J;

        X = dot(ELXY(:,1),SFL);
        Y = dot(ELXY(:,2),SFL);

        FX = FX0 + FX1*X;
        FY = FX0 + FY1*Y;
        %QZ = QZ0 + QZX1*X + QZY1*Y;
        QZ= DATA.QZ;
        QZ = QZ(X,Y);
        
        CNST = det(J)*FULL.GAUSWT(NG);

        ELK11 = ELK11 + A11*(GDSFL(:,1))*GDSFL(:,1)'*CNST +...
                        A66*(GDSFL(:,2))*GDSFL(:,2)'*CNST;
        ELK12 = ELK12 + A12*(GDSFL(:,1))*GDSFL(:,2)'*CNST +...
                        A66*(GDSFL(:,2))*GDSFL(:,1)'*CNST;
        
        ELK22 = ELK22 + A66*(GDSFL(:,1))*GDSFL(:,1)'*CNST +...
                        A22*(GDSFL(:,2))*GDSFL(:,2)'*CNST;
        %ELK33 = ELK33 + (A55*(GDSFL(:,1))*GDSFL(:,1)'+A44*(GDSFL(:,2))*GDSFL(:,2)'+...
        %                K*(SFL)*SFL')*CNST;
        ELK44 = ELK44 + (D11*GDSFL(:,1)*GDSFL(:,1)'+D66*GDSFL(:,2)*GDSFL(:,2)')*CNST;
        ELK45 = ELK45 + (D12*GDSFL(:,1)*GDSFL(:,2)'+D66*GDSFL(:,2)*GDSFL(:,1)')*CNST;
        ELK55 = ELK55 + (D66*GDSFL(:,1)*GDSFL(:,1)'+D22*GDSFL(:,2)*GDSFL(:,2)')*CNST;
                        

        ELF1  = ELF1 + F0*FX*SFL*CNST;
        ELF2 =  ELF2 + F0*FY*SFL*CNST;
        ELF3  = ELF3 + F0*QZ*SFL*CNST;
    end

    for NG = 1:size(REDU.GAUSWT,1)
        SFL   = REDU.SFL(:,NG);
        DSFL  = REDU.DSFL(:,:,NG);
        J = ELXY'*DSFL;
        GDSFL = DSFL/J;
        
        if(NONLIN)
            %U = dot(ELU,SFL);
            DWX = dot(ELU(3:NDF:end),GDSFL(:,1));
            DWY = dot(ELU(3:NDF:end),GDSFL(:,2));
            DUX = dot(ELU(1:NDF:end),GDSFL(:,1));
            DUY = dot(ELU(1:NDF:end),GDSFL(:,2));
            DVX = dot(ELU(2:NDF:end),GDSFL(:,1));
            DVY = dot(ELU(2:NDF:end),GDSFL(:,2));
           
        else
            U = 0;
            DWX = 0;
            DWY = 0;
        end
        CNST  = det(J)*REDU.GAUSWT(NG);
        
        ELK13 = ELK13 + 0.5*(GDSFL(:,1)*(A11*DWX*GDSFL(:,1)' + A12*DWY*GDSFL(:,2)')*CNST +...
                        A66*GDSFL(:,2)*(DWX*GDSFL(:,2)'+DWY*GDSFL(:,1)')*CNST);
        ELK23 = ELK23 + 0.5*(GDSFL(:,2)*(A12*DWX*GDSFL(:,1)' + A22*DWY*GDSFL(:,2)')*CNST +...
                        A66*GDSFL(:,1)*(DWX*GDSFL(:,2)'+DWY*GDSFL(:,1)')*CNST);
        ELK31 = ELK31 + GDSFL(:,1)*(A11*DWX*GDSFL(:,1)' + A66*DWY*GDSFL(:,2)')*CNST +...
                        GDSFL(:,2)*(A66*DWX*GDSFL(:,2)'+A12*DWY*GDSFL(:,1)')*CNST;
        ELK32 = ELK32 + GDSFL(:,1)*(A12*DWX*GDSFL(:,2)'+A66*DWY*GDSFL(:,1)')*CNST+...
                        GDSFL(:,2)*(A66*DWX*GDSFL(:,1)'+A22*DWY*GDSFL(:,2)')*CNST;
        ELK33 = ELK33 + (A55*(GDSFL(:,1))*GDSFL(:,1)'+A44*(GDSFL(:,2))*GDSFL(:,2)'+...
                        K*(SFL)*SFL')*CNST + 0.5*((A11*DWX^2+A66*DWY^2)*GDSFL(:,1)*GDSFL(:,1)'+...
                        (A66*DWX^2+A22*DWY^2)*GDSFL(:,2)*GDSFL(:,2)'+...
                        (A12+A66)*DWX*DWY*(GDSFL(:,1)*GDSFL(:,2)'+GDSFL(:,2)*GDSFL(:,1)'))*CNST;
        ELK34 = ELK34 + A55*GDSFL(:,1)*SFL'*CNST;
        ELK35 = ELK35 + A44*GDSFL(:,2)*SFL'*CNST;
        ELK44 = ELK44 + (A55*(SFL)*SFL')*CNST;
        ELK55 = ELK55 + A44*(SFL)*SFL'*CNST;
        
        PHIX    = dot(ELU(4:NDF:end),SFL);
        PHIY    = dot(ELU(5:NDF:end),SFL);

        DPHIXDX = dot(ELU(4:NDF:end),GDSFL(:,1));
        DPHIXDY = dot(ELU(4:NDF:end),GDSFL(:,2));

        DPHIYDX = dot(ELU(5:NDF:end),GDSFL(:,1));
        DPHIYDY = dot(ELU(5:NDF:end),GDSFL(:,2));

        [epsTop,sigTop,epsBot,sigBot,gamShear,tauShear] = POSTPRO_FSDT( ...
            DUX,DUY,DVX,DVY,DWX,DWY, ...
            PHIX,PHIY,DPHIXDX,DPHIXDY,DPHIYDX,DPHIYDY, ...
            DATA, H/2);

        STRS_TOP(NG,:) = [sigTop(:).' tauShear(:).'];
        STRS_BOT(NG,:) = [sigBot(:).' tauShear(:).'];

        %SGXX_BAR = (STRS_TOP)*10^2/7.8E6;

        if(NONLIN>1)

           Gx = GDSFL(:,1);
           Gy = GDSFL(:,2);  
           Nxx0 = A11*DUX + A12*DVY;
           Nyy0 = A12*DUX + A22*DVY;
           Nxy0 = A66*(DUY + DVX);
           
           TAN13 = TAN13 + 0.5*(GDSFL(:,1)*(A11*DWX*GDSFL(:,1)' + A12*DWY*GDSFL(:,2)')*CNST +...
                        A66*GDSFL(:,2)*(DWX*GDSFL(:,2)'+DWY*GDSFL(:,1)')*CNST);
           %TAN31 = TAN13';
           TAN23 = TAN23 + 0.5*(GDSFL(:,2)*(A12*DWX*GDSFL(:,1)' + A22*DWY*GDSFL(:,2)')*CNST +...
                        A66*GDSFL(:,1)*(DWX*GDSFL(:,2)'+DWY*GDSFL(:,1)')*CNST);
           %TAN32 = TAN23';
           % TAN33 = TAN33 + (A55*(GDSFL(:,1))*GDSFL(:,1)'+A44*(GDSFL(:,2))*GDSFL(:,2)' + ...
           %                  (A11*(DUX + 0.5*DWX^2)+A12*(DVY+0.5*DWY^2))*GDSFL(:,1)*GDSFL(:,1)' +...
           %                  (A12*(DUX + 0.5*DWX^2)+A22*(DVY+0.5*DWY^2))*GDSFL(:,2)*GDSFL(:,2)' +...
           %                  (A66*(DUY+DVX+DWX*DWY)*GDSFL(:,1)*(GDSFL(:,2)'))+...
           %                  (A12+A66)*DWX*DWY*((GDSFL(:,1)*GDSFL(:,2)')+(GDSFL(:,2)*GDSFL(:,1)')))*CNST-ELK33;
           TAN33 = TAN33 + ( ...
                            (Nxx0 + A11*DWX^2 + A66*DWY^2)*(Gx*Gx') + ...
                            (Nyy0 + A66*DWX^2 + A22*DWY^2)*(Gy*Gy') + ...
                            (Nxy0 + (A12 + A66)*DWX*DWY)*(Gx*Gy' + Gy*Gx') ...
                            )*CNST;
        end
    end
    ELK21 = ELK12';
    ELK43 = ELK34';
    ELK53 = ELK35';
    ELK54 = ELK45';

    ELK(1:NDF:end,1:NDF:end) = ELK11;
    ELK(1:NDF:end,2:NDF:end) = ELK12;
    ELK(1:NDF:end,3:NDF:end) = ELK13;
    ELK(1:NDF:end,4:NDF:end) = ELK14;
    ELK(1:NDF:end,5:NDF:end) = ELK15;
    ELK(2:NDF:end,1:NDF:end) = ELK21;
    ELK(2:NDF:end,2:NDF:end) = ELK22;
    ELK(2:NDF:end,3:NDF:end) = ELK23;
    ELK(2:NDF:end,4:NDF:end) = ELK24;
    ELK(2:NDF:end,5:NDF:end) = ELK25;
    ELK(3:NDF:end,1:NDF:end) = ELK31;
    ELK(3:NDF:end,2:NDF:end) = ELK32;
    ELK(3:NDF:end,3:NDF:end) = ELK33;
    ELK(3:NDF:end,4:NDF:end) = ELK34;
    ELK(3:NDF:end,5:NDF:end) = ELK35;
    ELK(4:NDF:end,1:NDF:end) = ELK41;
    ELK(4:NDF:end,2:NDF:end) = ELK42;
    ELK(4:NDF:end,3:NDF:end) = ELK43;
    ELK(4:NDF:end,4:NDF:end) = ELK44;
    ELK(4:NDF:end,5:NDF:end) = ELK45;
    ELK(5:NDF:end,1:NDF:end) = ELK51;
    ELK(5:NDF:end,2:NDF:end) = ELK52;
    ELK(5:NDF:end,3:NDF:end) = ELK53;
    ELK(5:NDF:end,4:NDF:end) = ELK54;
    ELK(5:NDF:end,5:NDF:end) = ELK55;
    ELF(1:NDF:end) = ELF1;
    ELF(2:NDF:end) = ELF2;
    ELF(3:NDF:end) = ELF3;
    if(NONLIN>1)
        TAN(1:NDF:end,3:NDF:end) = TAN13;
        %TAN(3:NDF:end,1:NDF:end) = TAN31;
        TAN(2:NDF:end,3:NDF:end) = TAN23;
        %TAN(3:NDF:end,2:NDF:end) = TAN32;
        TAN(3:NDF:end,3:NDF:end) = TAN33;
    end
    
    if(NONLIN>1)
        ELF = ELF - ELK*ELU;
        ELK = ELK + TAN;
    end
end

function [epsTop,sigTop,epsBot,sigBot,gamShear,tauShear] = POSTPRO_FSDT( ...
    DUX,DUY,DVX,DVY,DWX,DWY, ...
    PHIX,PHIY,DPHIXDX,DPHIXDY,DPHIYDX,DPHIYDY, ...
    DATA,z)

    Q11 = DATA.Q11;
    Q12 = DATA.Q12;
    Q21 = DATA.Q21;
    Q22 = DATA.Q22;
    Q44 = DATA.Q44;
    Q55 = DATA.Q55;
    Q66 = DATA.Q66;

    % membrane strains at mid-plane
    EXX0 = DUX + 0.5*DWX^2;
    EYY0 = DVY + 0.5*DWY^2;
    GXY0 = DUY + DVX + DWX*DWY;

    % curvatures
    KXX = DPHIXDX;
    KYY = DPHIYDY;
    KXY = DPHIXDY + DPHIYDX;

    % top and bottom strains
    EXX_T = EXX0 + z*KXX;
    EYY_T = EYY0 + z*KYY;
    GXY_T = GXY0 + z*KXY;

    EXX_B = EXX0 - z*KXX;
    EYY_B = EYY0 - z*KYY;
    GXY_B = GXY0 - z*KXY;

    % in-plane stresses
    sxx_top = Q11*EXX_T + Q12*EYY_T;
    syy_top = Q21*EXX_T + Q22*EYY_T;
    txy_top = Q66*GXY_T;

    sxx_bot = Q11*EXX_B + Q12*EYY_B;
    syy_bot = Q21*EXX_B + Q22*EYY_B;
    txy_bot = Q66*GXY_B;

    % transverse shear strains
    gxz = PHIX + DWX;
    gyz = PHIY + DWY;

    % transverse shear stresses
    txz = (DATA.A55 / DATA.H) * gxz;
    tyz = (DATA.A44 / DATA.H) * gyz;

    epsTop   = [EXX_T; EYY_T; GXY_T];
    sigTop   = [sxx_top; syy_top; txy_top];

    epsBot   = [EXX_B; EYY_B; GXY_B];
    sigBot   = [sxx_bot; syy_bot; txy_bot];

    gamShear = [gxz; gyz];
    tauShear = [txz; tyz];
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

