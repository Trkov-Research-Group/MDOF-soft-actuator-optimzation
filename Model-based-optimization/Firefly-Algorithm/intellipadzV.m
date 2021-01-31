                    

clc;
clear;
close all;

%% Problem Definition

CostFunctionz=@(model) mahsazV(model);


nVar=14;                 % Number of Decision Variables


VarMin_x=0;
VarMin_x2=(0.06339/2);
VarMin_L1= 0.06339;
VarMin_L11=0.006;
VarMin_b0=0.005;
VarMin_h=0.044;
VarMin_A11=0.006;
VarMin_b2=0.005;
VarMin_h2=0.0439999;
VarMin_A22=0.006;
VarMin_L2222=0.006;




% VarMin_x=0;
% VarMin_L1= 0.06339;
% VarMin_L11=0.001;
% VarMin_b0=0.0005;
% VarMin_h=0.003;
% VarMin_A11=0.0001;
% VarMin_b2=0.001;
% VarMin_h2=0.003;
% VarMin_A22=0.003;
% VarMin_L2222=0.0006;


% VarMin_L3=0.06339;
% VarMin_L33=0.01;
% VarMin_b3=0.004;
% VarMin_h3=0.044;
% VarMin_A33=0.02;
        


VarMax_x=0; 
VarMax_x2=0.06339/2;
VarMax_L1=0.06339;           
VarMax_L11=0.025;
VarMax_b0=0.025;
VarMax_h=0.044;
VarMax_A11=0.05;
VarMax_b2=0.025;
VarMax_h2=0.044;
VarMax_A22=0.05;
VarMax_L2222=0.025;


% VarMax_L3=0.06339;            
% VarMax_L33=0.02;
% VarMax_b3=0.01693333; 
% VarMax_h3=0.044;
% VarMax_A33=0.06;



%% Firefly Algorithm Parameters

MaxIt=500;         % Maximum Number of Iterations

nPop=25;            % Number of Fireflies (Swarm Size)

gamma=1;            % Light Absorption Coefficient

beta0=2;            % Attraction Coefficient Base Value

q=0.5;              % Intensification Factor (Selection Pressure)

alpha=0.2;          % Mutation Coefficient

alpha_damp=0.98;    % Mutation Coefficient Damping Ratio


delta_L1=0.05*(VarMax_L1-VarMin_L1);     % Uniform Mutation Range
delta_L11=0.05*(VarMax_L11-VarMin_L11);
delta_b0=0.05*(VarMax_b0-VarMin_b0);    % Uniform Mutation Range
delta_x=0.05*(VarMax_x-VarMin_x);
delta_x2=0.05*(VarMax_x2-VarMin_x2);
delta_h=0.05*(VarMax_h-VarMin_h);
delta_A11=0.05*(VarMax_A11-VarMin_A11);



delta_b2=0.05*(VarMax_b2-VarMin_b2);
delta_h2=0.05*(VarMax_h2-VarMin_h2);
delta_A22=0.05*(VarMax_A22-VarMin_A22);
delta_L2222=0.05*(VarMax_L2222-VarMin_L2222);

% delta_L3=0.05*(VarMax_L3-VarMin_L3);     % Uniform Mutation Range
% delta_L33=0.05*(VarMax_L33-VarMin_L33);
% delta_b3=0.05*(VarMax_b3-VarMin_b3);    % Uniform Mutation Range
% delta_h3=0.05*(VarMax_h3-VarMin_h3);
% delta_A33=0.05*(VarMax_A33-VarMin_A33);

     

m=2;

if isscalar(VarMin_L1) && isscalar(VarMax_L1)
    dmax_L1 = (VarMax_L1-VarMin_L1)*sqrt(nVar);
else
    dmax_L1 = norm(VarMax_L1-VarMin_L1);
end

if isscalar(VarMin_L11) && isscalar(VarMax_L11)
    dmax_L11 = (VarMax_L11-VarMin_L11)*sqrt(nVar);
else
    dmax_L11 = norm(VarMax_L11-VarMin_L11);
end

if isscalar(VarMin_x) && isscalar(VarMax_x)
    dmax_x = (VarMax_x-VarMin_x)*sqrt(nVar);
else
    dmax_x = norm(VarMax_x-VarMin_x);
end

if isscalar(VarMin_x2) && isscalar(VarMax_x2)
    dmax_x2 = (VarMax_x2-VarMin_x2)*sqrt(nVar);
else
    dmax_x2 = norm(VarMax_x2-VarMin_x2);
end


if isscalar(VarMin_h) && isscalar(VarMax_h)
    dmax_h = (VarMax_h-VarMin_h)*sqrt(nVar);
else
    dmax_h = norm(VarMax_h-VarMin_h);
end

if isscalar(VarMin_b0) && isscalar(VarMax_b0)
    dmax_b0 = (VarMax_b0-VarMin_b0)*sqrt(nVar);
else
    dmax_b0 = norm(VarMax_b0-VarMin_b0);
end

if isscalar(VarMin_A11) && isscalar(VarMax_A11)
    dmax_A11 = (VarMax_A11-VarMin_A11)*sqrt(nVar);
else
    dmax_A11 = norm(VarMax_A11-VarMin_A11);
end





if isscalar(VarMin_b2) && isscalar(VarMax_b2)
    dmax_b2 = (VarMax_b2-VarMin_b2)*sqrt(nVar);
else
    dmax_b2 = norm(VarMax_b2-VarMin_b2);
end

if isscalar(VarMin_h2) && isscalar(VarMax_h2)
    dmax_h2 = (VarMax_h2-VarMin_h2)*sqrt(nVar);
else
    dmax_h2 = norm(VarMax_h2-VarMin_h2);
end


if isscalar(VarMin_A22) && isscalar(VarMax_A22)
    dmax_A22 = (VarMax_A22-VarMin_A22)*sqrt(nVar);
else
    dmax_A22 = norm(VarMax_A22-VarMin_A22);
end

if isscalar(VarMin_L2222) && isscalar(VarMax_L2222)
    dmax_L2222 = (VarMax_L2222-VarMin_L2222)*sqrt(nVar);
else
    dmax_L2222 = norm(VarMax_L2222-VarMin_L2222);
end




%% Initialization

% Empty Firefly Structure
fireflyz.Position=[];
fireflyz.Cost=[];

popz=repmat(fireflyz,nPop,1);


% Initialize Best Solution Ever Found
BestSolz.Cost= inf;




% Create Initial Fireflies
for i=1:nPop
    model=CreateModelYo11V(); 
    popz(i).Position=model;
    popz(i).Cost=CostFunctionz(popz(i).Position);
    
    if popz(i).Cost<=BestSolz.Cost
        BestSolz=popz(i);
    end
end
BestCostz=zeros(MaxIt,1);
BestCostzYo2=zeros(MaxIt,1);
BestCostzYo3=zeros(MaxIt,1);
BestCostz1=zeros(MaxIt,1);
BestCostz2=zeros(MaxIt,1);
BestCostz3=zeros(MaxIt,1);
BestCostz4=zeros(MaxIt,1);
BestCostz5=zeros(MaxIt,1);
BestCostz6=zeros(MaxIt,1);
BestCostz7=zeros(MaxIt,1);
BestCostz8=zeros(MaxIt,1);
BestCostz9=zeros(MaxIt,1);
BestCostz10=zeros(MaxIt,1);
BestCostz11=zeros(MaxIt,1);
BestCostz12=zeros(MaxIt,1);
BestCostz_delta=zeros(MaxIt,1);


rij.L1=[];
rij.L11=[];
rij.x=[];
rij.x2=[];
rij.h=[];
rij.b0=[];
rij.A11=[];
rij.b2=[];
rij.h2=[];
rij.A22=[];
rij.L2222=[];
rij.L3=[];
rij.L33=[];
rij.h3=[];
rij.L333=[];
rij.b3=[];
rij.A33=[];

repmat(rij,nPop,1);

beta.L1=[];
beta.L11=[];
beta.x=[];
rij.x2=[];
beta.h=[];
beta.b0=[];
beta.A11=[];
beta.b2=[];
beta.h2=[];
beta.A22=[];
beta.L2222=[];
beta.L3=[];
beta.L33=[];
beta.h3=[];
beta.L333=[];
beta.b3=[];
beta.A33=[];

repmat(beta,nPop,1);



for it=1:MaxIt
    
    newpopz=repmat(fireflyz,nPop,1);
  

    for i=1:nPop
        newpopz(i).Cost = inf; 
        for j=1:nPop    
            if popz(j).Cost <= popz(i).Cost
                
                
                
                rij.A11=norm(popz(i).Position.A11-popz(j).Position.A11)/dmax_A11;
                rij.L1=norm(popz(i).Position.L1-popz(j).Position.L1)/dmax_L1;
                rij.L11=norm(popz(i).Position.L11-popz(j).Position.L11)/dmax_L11;
                rij.x=norm(popz(i).Position.x-popz(j).Position.x)/dmax_x;
                rij.x2=norm(popz(i).Position.x2-popz(j).Position.x2)/dmax_x2;
                rij.h=norm(popz(i).Position.h-popz(j).Position.h)/dmax_h;
                rij.b0=norm(popz(i).Position.b0-popz(j).Position.b0)/dmax_b0;
                rij.b2=norm(popz(i).Position.b2-popz(j).Position.b2)/dmax_b2;
                rij.h2=norm(popz(i).Position.h2-popz(j).Position.h2)/dmax_h2;
                rij.A22=norm(popz(i).Position.A22-popz(j).Position.A22)/dmax_A22;
                rij.L2222=norm(popz(i).Position.L2222-popz(j).Position.L2222)/dmax_L2222;
                
                beta.L1=beta0*exp(-gamma*rij.L1^m);
                beta.L11=beta0*exp(-gamma*rij.L11^m);
                beta.x=beta0*exp(-gamma*rij.x^m);
                beta.x2=beta0*exp(-gamma*rij.x2^m);
                beta.h=beta0*exp(-gamma*rij.h^m);
                beta.b0=beta0*exp(-gamma*rij.b0^m);
                beta.A11=beta0*exp(-gamma*rij.A11^m);
                beta.b2=beta0*exp(-gamma*rij.b2^m);
                beta.h2=beta0*exp(-gamma*rij.h2^m);
                beta.A22=beta0*exp(-gamma*rij.A22^m);
                beta.L2222=beta0*exp(-gamma*rij.L2222^m);
                
                
                e_L1=delta_L1*unifrnd(0,+1);
                e_L11=delta_L11*unifrnd(0,+1);
                e_x=delta_x*unifrnd(0,+1);
                e_x2=delta_x2*unifrnd(0,+1);
                e_h=delta_h*unifrnd(0,+1);
                e_b0=delta_b0*unifrnd(0,+1);
                e_A11=delta_A11*unifrnd(0,+1);
                e_b2=delta_b2*unifrnd(0,+1);
                e_h2=delta_h2*unifrnd(0,+1);
                e_A22=delta_A22*unifrnd(0,+1);
                e_L2222=delta_L2222*unifrnd(0,+1);
                
                
                newsolz.Position.L1 = popz(i).Position.L1 ...
                    + beta.L1*rand.*(popz(j).Position.L1-popz(i).Position.L1) ...
                    + alpha*e_L1;
                newsolz.Position.L11 = popz(i).Position.L11 ...
                    + beta.L11*rand.*(popz(j).Position.L11-popz(i).Position.L11) ...
                    + alpha*e_L11;
                newsolz.Position.x = popz(i).Position.x ...
                    + beta.x*rand.*(popz(j).Position.x-popz(i).Position.x) ...
                    + alpha*e_x;
                newsolz.Position.x2 = popz(i).Position.x2 ...
                    + beta.x*rand.*(popz(j).Position.x-popz(i).Position.x) ...
                    + alpha*e_x;
                newsolz.Position.h = popz(i).Position.h ...
                    + beta.h*rand.*(popz(j).Position.h-popz(i).Position.h) ...
                    + alpha*e_h;
                newsolz.Position.b0 = popz(i).Position.b0 ...
                    + beta.b0*rand.*(popz(j).Position.b0-popz(i).Position.b0) ...
                    + alpha*e_b0;
                newsolz.Position.A11 = popz(i).Position.A11 ...
                    + beta.A11*rand.*(popz(j).Position.A11-popz(i).Position.A11) ...
                    + alpha*e_A11;
                
                newsolz.Position.L2 = 0.06339;
                
                newsolz.Position.b2 = popz(i).Position.b2 ...
                    + beta.b2*rand.*(popz(j).Position.b2-popz(i).Position.b2) ...
                    + alpha*e_b2;
                newsolz.Position.h2 = popz(i).Position.h2 ...
                    + beta.h2*rand.*(popz(j).Position.h2-popz(i).Position.h2) ...
                    + alpha*e_h2;
                newsolz.Position.A22 = popz(i).Position.A22 ...
                    + beta.A22*rand.*(popz(j).Position.A22-popz(i).Position.A22) ...
                    + alpha*e_A22;
                
                newsolz.Position.L2222 = popz(i).Position.L2222 ...
                    + beta.L2222*rand.*(popz(j).Position.L2222-popz(i).Position.L2222) ...
                    + alpha*e_L2222;
                

%                 
%                 
                if (newsolz.Position.b2) < (newsolz.Position.L2222)
                    minLb = newsolz.Position.L2222;
                    newsolz.Position.L2222 = newsolz.Position.b2;
                    newsolz.Position.b2 = minLb;
                end

% 
%                 if (newsolz.Position.L11) < (newsolz.Position.b0)
%                     minLb = newsolz.Position.L11;
%                     newsolz.Position.L11 = newsolz.Position.b0;
%                     newsolz.Position.b0 = minLb;
%                 end



                

                newsolz.Position.L1=max(newsolz.Position.L1,VarMin_L1);
                newsolz.Position.L1=min(newsolz.Position.L1,VarMax_L1);

                newsolz.Position.L11=max(newsolz.Position.L11,VarMin_L11);
                newsolz.Position.L11=min(newsolz.Position.L11,VarMax_L11);

                
                newsolz.Position.x=max(newsolz.Position.x,VarMin_x);
                newsolz.Position.x=min(newsolz.Position.x,VarMax_x);          

                newsolz.Position.x2=max(newsolz.Position.x2,VarMin_x2);
                newsolz.Position.x2=min(newsolz.Position.x2,VarMax_x2);
                
                newsolz.Position.h=max(newsolz.Position.h,VarMin_h);
                newsolz.Position.h=min(newsolz.Position.h,VarMax_h);

                
                newsolz.Position.b0=max(newsolz.Position.b0,VarMin_b0);
                newsolz.Position.b0=min(newsolz.Position.b0,VarMax_b0);

                
                newsolz.Position.A11=max(newsolz.Position.A11,VarMin_A11);
                newsolz.Position.A11=min(newsolz.Position.A11,VarMax_A11);
                newsolz.Position.L111=(((newsolz.Position.L1^2)+((newsolz.Position.L11-newsolz.Position.b0)^2)).^0.5);
                
                
                
                
                
                
                
                newsolz.Position.L2=0.06339;
                newsolz.Position.b2=max(newsolz.Position.b2,VarMin_b2);
                newsolz.Position.b2=min(newsolz.Position.b2,VarMax_b2);
                
                newsolz.Position.h2=max(newsolz.Position.h2,VarMin_h2);
                newsolz.Position.h2=min(newsolz.Position.h2,VarMax_h2);
                
                newsolz.Position.A22=max(newsolz.Position.A22,VarMin_A22);
                newsolz.Position.A22=min(newsolz.Position.A22,VarMax_A22);
                
                newsolz.Position.L2222=max(newsolz.Position.L2222,VarMin_L2222);
                newsolz.Position.L2222=min(newsolz.Position.L2222,VarMax_L2222);
                
                newsolz.Position.L222=((newsolz.Position.L2^2)+((((newsolz.Position.b2)/2)-((newsolz.Position.L2222)/2))^2))^(0.5);
                
                


                newsolz.Position.L3=newsolz.Position.L1;
                newsolz.Position.L33=newsolz.Position.L11;
                newsolz.Position.h3=newsolz.Position.h;
                newsolz.Position.b3=newsolz.Position.b0;
                newsolz.Position.A33=newsolz.Position.A11;
                newsolz.Position.L333=newsolz.Position.L111;
                newsolz.Position.d1=(0.0508-(2*newsolz.Position.b0)-newsolz.Position.b2);
                newsolz.Position.d11=(0.0508-(2*newsolz.Position.L11)-newsolz.Position.L2222);
%                 newsolz.Position.Ks=100000000000000;
                newsolz.Position.Ks=(2*6.72087912)/(0.0508-newsolz.Position.b0-newsolz.Position.b2-newsolz.Position.b3);
                
                
                L11=newsolz.Position.L11;
                x=newsolz.Position.x; 
                x2=newsolz.Position.x2;  
                h=newsolz.Position.h;
                L1=newsolz.Position.h;
                b0=newsolz.Position.b0;
                A11=newsolz.Position.A11;
                L111=newsolz.Position.L111;
                
                
                L2=newsolz.Position.L2;
                b2=newsolz.Position.b2;
                h2=newsolz.Position.h2;
                A22=newsolz.Position.A22;
                L2222=newsolz.Position.L2222;
                L222=newsolz.Position.L222;
                
                L3=newsolz.Position.L3;
                L33=newsolz.Position.L33;
                h3=newsolz.Position.h3;
                b3=newsolz.Position.b3;
                A33=newsolz.Position.A33;
                L333=newsolz.Position.L333;
                
                d1=newsolz.Position.d1;
                d11=newsolz.Position.d11;
                Ks=newsolz.Position.Ks;
                
                cos1=(L1/L111);
                cos3=(L3/L333);
                
                
                p1=-(103421);
                p3=(103421);
                p2=(-p1+p3);
                
                
                
                
                A1=(h*A11);
                q1=((p1*(A1*cos1))/L1);
                alpha=((L11-b0)/L1);
                E=885523.3145;
                R=(q1/((alpha^4)*E*h));
                F1=q1*L1;
                
                C2=(3*R*b0);
                C4=((3/4)*R*(b0*b0));
                C3=(C4/(L11*L11))-(C2/L11)-((1.5)*R*((log(L11))+1));
                C1=-((C2*(log(L11)))+(C3*L11)+(C4/L11)+(1.5*R*(log(L11))*L11));
                
                Y1=(C1+ (C2*(log((alpha*x)+b0))) + (C3*((alpha*x)+b0)) + (C4/((alpha*x)+b0)) + (1.5*R*(log((alpha*x)+b0))*((alpha*x)+b0)));
                Y11=(C1+ (C2*(log((alpha*x2)+b0))) + (C3*((alpha*x2)+b0)) + (C4/((alpha*x2)+b0)) + (1.5*R*(log((alpha*x2)+b0))*((alpha*x2)+b0)));
                
                
                
                cos2=(L2/L222);
                A2=(h2*A22);
                q2=((p2*(A2*cos2))/L2);
                alpha2=((L2222-b2)/L2);
                
                
                R2=(q2/((alpha2^4)*E*h2));
                F2=q2*L2;
                
                C22=((6/4)*R2*b2);
                C42=((3/8)*R2*(b2*b2));
                C32=(C42/(L2222*L2222))-(C22/L2222)-((3/4)*R2*((log(L2222))+1));
                C12=-((C22*(log(L2222)))+(C32*L2222)+(C42/L2222)+((3/4)*R2*(log(L2222))*L2222));
                
                Y2=(C12+ (C22*(log((alpha2*x)+b2))) + (C32*((alpha2*x)+b2)) + (C42/((alpha2*x)+b2)) + ((3/4)*R2*(log((alpha2*x)+b2))*((alpha2*x)+b2)));
                Y22=(C12+ (C22*(log((alpha2*x2)+b2))) + (C32*((alpha2*x2)+b2)) + (C42/((alpha2*x2)+b2)) + ((3/4)*R2*(log((alpha2*x2)+b2))*((alpha2*x2)+b2)));
                
                A3=(h*A11);
                q3=((p3*(A3*cos1))/L1);
                alpha3=((L11-b0)/L1);
                R3=(q3/((alpha^4)*E*h));
                F3=q3*L1;
                
                C23=(3*R3*b3);
                C43=((3/4)*R3*(b3*b3));
                C33=(C43/(L33*L33))-(C23/L33)-((1.5)*R3*((log(L33))+1));
                C13=-((C23*(log(L33)))+(C33*L33)+(C43/L33)+(1.5*R3*(log(L33))*L33));
                
                Y3=-(C13+ (C23*(log((alpha3*x)+b3))) + (C33*((alpha3*x)+b3)) + (C43/((alpha3*x)+b3)) + (1.5*R3*(log((alpha3*x)+b0))*((alpha3*x)+b3)));
                Y33=-(C13+ (C23*(log((alpha3*x2)+b3))) + (C33*((alpha3*x2)+b3)) + (C43/((alpha3*x2)+b3)) + (1.5*R3*(log((alpha3*x2)+b0))*((alpha3*x2)+b3)));
                
                
                K1=abs(F1/Y1);
                K2=abs(F2/Y2);
                K3=abs(F3/Y3);
                K11=abs(F1/Y11);
                K22=abs(F2/Y22);
                K33=abs(F3/Y33);
                
                
                Yo2=((Y2)+((Ks*Y1)/(K2*(1+(Ks/K1))))-((Ks*K3*Y3)/(K2*(K3+Ks))))/(1-((Ks*Ks)/(K2*(K3+Ks)))-((Ks*Ks)/(K2*(K1+Ks)))+((2*Ks)/K2));
                Yo22=((Y22)+((Ks*Y11)/(K22*(1+(Ks/K11))))-((Ks*K33*Y33)/(K22*(K33+Ks))))/(1-((Ks*Ks)/(K2*(K33+Ks)))-((Ks*Ks)/(K22*(K11+Ks)))+((2*Ks)/K22));
                
                Yo1=((((Ks/K1)*(Yo2))+(Y1))/(1+(Ks/K1)));
                Yo11=((((Ks/K11)*(Yo22))+(Y11))/(1+(Ks/K11)));
                
                Yo3=((((Ks/K3)*(Yo2))-(Y3))/(1+(Ks/K3))); 
                Yo33=((((Ks/K33)*(Yo22))-(Y33))/(1+(Ks/K33)));
                
                
                miu = 1;
                miu1 = 1;
                miu2 = 1;
                miu3 = 1;
                
                delta_L = miu2*((19994.8*0.0508*0.0508*L1)/((2*b0+L2222)*E*h2)).^2;
                
                max1= miu1*((max(0,-((2*b0)+L2222-0.03))).^2);
                max_delta = miu2*((((delta_L))).^2);
                max2=miu3*((max(0,(0.0508-(2*b0)+b2+d1))).^2);
                max3=miu*((max(0,-(0.0508-(2*b0)+b2+d1))).^2);
                max4=miu*((max(0,(-d1))).^2);
                max5=miu*((max(0,(-d11))).^2);
                max6=miu*((max(0,(0.0508-(2*L11)+L2222+d11))).^2);
                max7=miu*((max(0,(-(0.0508-(2*L11)+L2222+d11)))).^2);
                max8=miu*((max(0,(-((b0/L11)-0.1)))).^2);
                max9=miu*((max(0,((b0/L11)-0.60))).^2);
                max10=miu*((max(0,(-((L2222/b2)-0.1)))).^2);
                max11=miu*((max(0,((L2222/b2)-0.60))).^2);
                
                newsolz.Position.Yo2=Yo2;
                newsolz.Position.Yo1=Yo1;
                newsolz.Position.Yo3=Yo3;
                newsolz.Position.Y1=Y1;
                newsolz.Position.Y2=Y2;
                newsolz.Position.Y3=Y3;
                newsolz.Position.Yo22=Yo22;
                newsolz.Position.Yo11=Yo11;
                newsolz.Position.Yo33=Yo33;
                newsolz.Position.Y11=Y11;
                newsolz.Position.Y22=Y22;
                newsolz.Position.Y33=Y33;
                                
                newsolz.Cost=CostFunctionz(newsolz.Position);

                if newsolz.Cost <= newpopz(i).Cost 
                    
                    newpopz(i) = newsolz;
                    
                    if newpopz(i).Cost<=BestSolz.Cost
                        BestSolz=newpopz(i);
                        
                    end
                end
                
            end
        end
        
    end
    
    
    
    % Merge
    popz=[popz
        newpopz];  %#ok
    
    % Sort
    [~, SortOrder]=sort([popz.Cost]);
    popz=popz(SortOrder);
    
    % Truncate
    popz=popz(1:nPop);
    
    % Store Best Cost Ever Found
    BestCostz(it)=BestSolz.Cost;
    BestCostzYo2(it)=Yo2;
    BestCostzYo3(it)=Yo3;
    BestCostz1(it)=Y1;
    BestCostz2(it)=Y2;
    BestCostz3(it)=Y3;
    BestCostz5(it)=max5;
    BestCostz6(it)=max6;
    BestCostz7(it)=max7;
    BestCostz8(it)=max8;
    BestCostz9(it)=max9;
    BestCostz10(it)=max10;
    BestCostz11(it)=max11;
    BestCostz_delta(it)=max_delta;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost z = ' num2str(BestCostz(it))]);
    
    % Damp Mutation Coefficient
    alpha = alpha*alpha_damp;
end

%% Results

figure;
%plot(BestCost,'LineWidth',2);
semilogy(-BestCostz,'LineWidth',2);
xlabel('Iteration');
ylabel('Cost');
xlim([-1,50]);
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostzYo2,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost Yo2');
% grid on;
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostzYo3,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost Yo3');
% grid on;
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz1,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost 1');
% grid on;
% 
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz2,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 2');
% grid on;


   
% % L1 = newsolz.Position.L1;   % height 
% % b0 = newsolz.Position.b0 ;  % top side
% % L11 = newsolz.Position.L11 ;   % base 
% % %Frame vertices
% % A = [0 0] ;
% % B = [L11 0] ;
% % C = [b0 L1] ;
% % D = [0 L1] ;  
% % coor = [A ; B; C; D] ;  
% % patch(coor(:,1), coor(:,2),'r')
% % 
% % L2 = newsolz.Position.L2;   % height 
% % b2 = newsolz.Position.b2 ;  % top side
% % L2222 = newsolz.Position.L2222 ;   % base 
% % %%Frame vertices
% % A = [0+L11+d11/2 0] ;
% % B = [L2222+L11+d11/2 0] ;
% % % C = [x y] ;
% % C = [(0.5*(L2222-b2)+b2)+L11+d11/2 L2] ;
% % D = [0.5*(L2222-b2)+L11+d11/2 L2] ;  
% % coor = [A ; B; C; D] ;  
% % patch(coor(:,1), coor(:,2),'r')
% % 
% % 
% % L3 = newsolz.Position.L3;   % height 
% % b3 = newsolz.Position.b3 ;  % top side
% % L33 = newsolz.Position.L33 ;   % base 
% % % L1 = 10 ;  % height 
% % % b0 = 4 ;  % top side
% % % L11 = 8 ;   % base 
% % %Frame vertices
% % A = [0+L11+d11+L2222 0] ;
% % B = [L33+L11+d11+L2222 0] ;
% % C = [L33+L11+d11+L2222 L3] ;
% % D = [L33+L11+d11+L2222-b3 L3] ;  
% % coor = [A ; B; C; D] ;  
% % patch(coor(:,1), coor(:,2),'r')


% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz3,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 3');
% grid on;

% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz5,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 5');
% grid on;
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz6,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 6');
% grid on;
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz7,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 7');
% grid on;
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz8,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 8');
% grid on;
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz9,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 9');
% grid on;
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz10,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 10');
% grid on;
% 
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz11,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z 11');
% grid on;
% 
% 
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCostz_delta,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost z _ delta');
% grid on;