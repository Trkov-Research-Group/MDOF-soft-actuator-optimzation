function model=CreateModelYo11V()


      x=unifrnd(0,0.06339/2);  
      x2=unifrnd(0.06339/2,0.04);
      L11=unifrnd(0.004,0.0169);
      L1=0.06339;    
      b0=unifrnd(0.004,0.0169);
      h=unifrnd(0.004,0.044);
      A11=unifrnd(0.004,0.05); 
      L111=(((L1*L1)+((L11-b0)*(L11-b0)))^0.5);

      L2=0.06339;
      b2=unifrnd(0.004,0.0169);
      h2=unifrnd(0.004,0.044);
      A22=unifrnd(0.004,0.05);
      L2222=unifrnd(0.004,0.0169);
      L222=(((L2*L2)+((b2-L2222)*(b2-L2222)))^0.5);

      b3=b0;

      Ks=(2*6.72087912)/(0.0508-b0-b2-b3);
      

      model.L1=L1;
      model.L11=L11;
      model.L111=L111;
      model.b0=b0;
      model.h=h;
      model.x=x;
      model.x2=x2;
      model.A11=A11;
      
      
      model.L3=L1;
      model.L33=L11;
      model.L333=L111;
      model.b3=b0;
      model.h3=h;
      model.A33=A11;
      
      L3=model.L3;
      L33=model.L33;
      L333=model.L333;
      b3=model.b3;
      h3=model.h3;
      A33=model.A33;

      model.L2=L2;
      model.b2=b2;
      model.h2=h2;
      model.A22=A22;
      model.L2222=L2222;
      

      model.Ks=Ks;
      model.d1=(0.0254-model.b0-(model.b2/2));
      model.d11=(0.0254-model.L11-(model.L2222/2));
      
      
end