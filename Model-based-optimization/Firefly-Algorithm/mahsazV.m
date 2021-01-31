function z = mahsazV(model)
     

      
        x=model.x;
        x2=model.x2;
        L1=model.L1;
        L11=model.L11;
        L111=model.L111;
        b0=model.b0;
        h=model.h;
        A11=model.A11;
        
        L2=model.L2;
        L2222=model.L2222;
        b2=model.b2;
        h2=model.h2;
        A22=model.A22;
        L222=((L2*L2)+((b2-L2222)*(b2-L2222)))^(0.5);
        
        L33=model.L33;
        b3=model.b3;
        d1=(0.0548-((2*b0)+b2));
        d11=(0.0548-((2*L11)+L2222));
        Ks=model.Ks;     
        cos1=(L1/L111);


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
        
        
        miu = 1;
        miu1 = 1;
        miu2 = 1;
        miu3 = 1;
        
        Yo2=((Y2)+((Ks*Y1)/(K2*(1+(Ks/K1))))-((Ks*K3*Y3)/(K2*(K3+Ks))))/(1-((Ks*Ks)/(K2*(K3+Ks)))-((Ks*Ks)/(K2*(K1+Ks)))+((2*Ks)/K2));
        Yo22=((Y22)+((Ks*Y11)/(K22*(1+(Ks/K11))))-((Ks*K33*Y33)/(K22*(K33+Ks))))/(1-((Ks*Ks)/(K2*(K33+Ks)))-((Ks*Ks)/(K22*(K11+Ks)))+((2*Ks)/K22));
        Yo1=-((((Ks/K1)*(Yo2))+(Y1))/(1+(Ks/K1)));
        Yo33=((((Ks/K33)*(Yo22))-(Y33))/(1+(Ks/K33)));
  
        
        
        z=  -(Yo1) + ...
            -miu1*10*(abs(Yo33)) + ...
            -miu1*1*(Y33) + ...
            -miu1*1*(Y11) + ...
            miu1*1*((max(0,(((((((Ks/K1)*(Yo2))+(Y1))/(1+(Ks/K1)))-b0)+(Yo2+b2/2))-d1))).^2) + ...
            miu1*10*((max(0,-((2*b0)+L2222-0.03))).^2) + ...
            miu3*10*((max(0,(0.0508-((2*abs(b0))+abs(b2)+abs(d1))))).^2) + ...
            miu3*10*((max(0,-(0.0508-((2*abs(b0))+abs(b2)+abs(d1))))).^2) + ...
            miu3*10*((max(0,(0.0508-((2*abs(L11))+abs(L2222)+abs(d11))))).^2) + ...
            miu*10*((max(0,(-(0.0508-((2*abs(L11))+abs(L2222)+abs(d11)))))).^2) + ...
            miu2*100*((max(0,(-d1))).^2) + ...
            miu2*100*((max(0,(-d11))).^2) + ...
            miu*((max(0,(-((b0/L11)-0.1)))).^2) + ...
            miu*100*((max(0,((b0/L11)-0.9))).^2) + ...
            miu*((max(0,(-((L2222/b2)-0.1)))).^2) + ...
            miu*100*((max(0,((L2222/b2)-0.4))).^2);

        
        
end