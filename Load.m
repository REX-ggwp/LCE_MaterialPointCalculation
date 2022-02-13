clear all
tic

dt=0.01; %Time Step
rate=0.1;%Applied Strain Rate
%dt=.01; 
%Material Parameters
MaterialProperties.Q=0.5; %Order parameter
MaterialProperties.JmEQ=5.7; %Locking stretch for equilibrium branch
MaterialProperties.JmNEQ=5.7;%Loacking stretch for nonequilibriun branch, can be an array
MaterialProperties.muEQ  = 0.25; %Shear modulus for equilibrium branch
MaterialProperties.muNEQ = 1.25; %Shear modulus for nonequilibrium branch, can be an array
MaterialProperties.eta_N = 800;  %Viscosity, can be an array
MaterialProperties.type = {'NG'};%Type of spring for nonequilibriun branch, can be an array,'NG'=Neo Gent, 'NC'=Neo Classical
MaterialProperties.eta_D = 40; %Director Viscosity
MaterialProperties.ks = 1/0.7; %Stress factor for the Director Viscosity

%{
load('Spectrum24.mat')
MaterialProperties.muEQ=Spectrum.mueq; 
MaterialProperties.muNEQ= Spectrum.muneq;
MaterialProperties.eta_N = Spectrum.eta*2*0.01;
MaterialProperties.type = {'NG','NG','NG','NG','NG','NG',...
                            'NG','NG','NG','NG','NG','NG',...
                            'NG','NG','NG','NG','NG','NG',...
                           'NG','NG','NG','NG','NG','NG'};
%}
%{
MaterialProperties.muEQ=0.432094506901689; 
MaterialProperties.muNEQ= [0.120882014698349, 0.149356619973996,  0.158309399513377, 0.207875052709189,...
                           0.245335300827003, 0.319215245328831,  0.657240162391083, 1.791235360992085...
                           5.427310645998254, 16.394293062514066, 42.034433575228036,27.315437165422267];
MaterialProperties.eta_N = [1.208820146983487e+05,  2.798670367268106e4, 5.558546550199304e+03, 1.367679037548190e+03, ...
                            3.024608256946087e+02,  73.742861895346550,  28.450357673312627,    14.529262265013752,...
                            8.249029601481190,      4.669153467968824,   2.243251293451123,     0.273154371654223]*2*0.01;
MaterialProperties.type = {'NG','NG','NG','NG','NG','NG',...
                           'NG','NG','NG','NG','NG','NG'};
%}
%{
MaterialProperties.muEQ=0.1054; 
MaterialProperties.muNEQ= [0.0047,0.0150,0.0405,0.0938,0.1783,0.3006,0.5131,1.01783,...
                           2.4372,6.5587,17.6547,40.5306,71.0178,85.3125,66.4163,78.1615];
MaterialProperties.eta_N = [47.4429,32.3189,18.8021,9.3796,3.8420,1.3952,0.5131,0.2193,...
                            0.1131,0.0656,0.0380,0.0188,0.0071,0.0018,3.0828E-04,7.8161E-05];
MaterialProperties.type = {'NG','NG','NG','NG','NG','NG','NG','NG',...
                            'NG','NG','NG','NG','NG','NG','NG','NG'};
%MaterialProperties.type = {'NC','NC','NC','NC','NC','NC','NC','NC',...
%                            'NC','NC','NC','NC','NC','NC','NC','NC'};
%MaterialProperties.type = {'NG'};
%}



theta0 = 90/180*pi; %Initial director angle

%Construct the object
LCEobj = LCE(MaterialProperties, theta0); %Initialize the object


%UniaxialTension
MS = 5.25;  %Maximum stretch
l1 = (1+dt*rate) : (dt*rate) : MS; %Loading, constant engineering strain rate
l1 =[l1,MS-dt*rate:-dt*rate:1];%Unloading

%Result is the structure to be saved
Result.MaterialProperties=MaterialProperties;
Result.LoadingRate=rate;
Result.theta0 = theta0;
%Loop over each step
for ii=1:length(l1)
    %ii th Step Starts, Update History Variable
    LCEobj.UpdateHistory;
    %Make Initial Guesses
    InitialGuess = LCEobj.CurrentState.F;
    InitialGuess = [InitialGuess(2,2),InitialGuess(1,2)+0.2*rate*dt*rand()];
    %Solve Boundary Conditions
    fun = @(x)ResStress(l1(ii),x, dt, LCEobj);
    [Solution,~,flag,~] = fsolve(fun, InitialGuess );%Solve for residue=0
    if flag<1
        error("fsolve error when solve F");
    end
    %Save all the variables into the structure Result
    Result.S11(ii) = LCEobj.S(1,1);
    Result.S22(ii) = LCEobj.S(2,2);
    Result.S21(ii) = LCEobj.S(2,1);
    Result.S12(ii) = LCEobj.S(1,2);
    Result.S11EQ(ii)  = LCEobj.S_eq(1,1);
    Result.S11NEQ(ii) = LCEobj.S(1,1)-LCEobj.S_eq(1,1);
    Result.p(ii)   = LCEobj.p;
    Result.theta(ii) = LCEobj.CurrentState.theta;
    Result.l1(ii)  = LCEobj.CurrentState.F(1,1);
    Result.l2(ii)  = LCEobj.CurrentState.F(2,2);
    Result.l3(ii)  = LCEobj.CurrentState.l3;
    Result.gamma(ii,:)  = LCEobj.CurrentState.F(1,2,:);
    Result.l1v(ii,:)    = reshape(LCEobj.CurrentState.Fv(1,1,:),[1,LCEobj.num_Spec]);
    Result.l2v(ii,:)    = reshape(LCEobj.CurrentState.Fv(2,2,:),[1,LCEobj.num_Spec]);
    Result.l3v(ii,:)    = LCEobj.CurrentState.l3v;
    Result.Fv12(ii,:)   = reshape(LCEobj.CurrentState.Fv(1,2,:),[1,LCEobj.num_Spec]);
    Result.Fv21(ii,:)   = reshape(LCEobj.CurrentState.Fv(2,1,:),[1,LCEobj.num_Spec]);
    Result.l1e(ii,:)    = reshape(LCEobj.Fe(1,1,:),[1,LCEobj.num_Spec]);
    Result.l2e(ii,:)    = reshape(LCEobj.Fe(2,2,:),[1,LCEobj.num_Spec]);
    Result.l3e(ii,:)    = LCEobj.l3e;
    Result.Fe12(ii,:)   = reshape(LCEobj.l3e,[1,LCEobj.num_Spec]);
    Result.Fe21(ii,:)   = reshape(LCEobj.l3e,[1,LCEobj.num_Spec]);
    Result.D_Dir(ii)    = LCEobj.D_Dir;
    Result.D_Net(ii)    = LCEobj.D_Net;
    Result.Psi(ii)      = LCEobj.Psi;
end
%clear dt flag fun ii InitialGuess l1 LCEobj MaterialProperties MS rate Solution theta0
toc

function Res = ResStress(l1, x, dt, obj)
    F = [l1,x(2);
         0, x(1);];
    obj.ApplyLoad(F,dt);
    s = obj.S;
    Res = [s(2,2)/s(1,1);s(1,2)/s(1,1)]; %Return desired residue
end

