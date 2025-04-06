function 	SE=Ks(n,NOCs,XofNs,Material_Elasstic,materils)

num_mater= materils(n);
D=cell2mat(Material_Elasstic(num_mater,6));
 
      XNI(1, 1) = -0.57735026919;
	  XNI(2, 1) = -0.57735026919;
      XNI(1, 2) =  0.57735026919;
	  XNI(2, 2) = -0.57735026919;
      XNI(1, 3) =  0.57735026919;
	  XNI(2, 3) =  0.57735026919;
      XNI(1, 4) = -0.57735026919;
	  XNI(2, 4) =  0.57735026919;

      X1 = XofNs(NOCs(n, 1), 1);
      Y1 = XofNs(NOCs(n, 1), 2);
      X2 = XofNs(NOCs(n, 2), 1);
      Y2 = XofNs(NOCs(n, 2), 2);
      X3 = XofNs(NOCs(n, 3), 1);
      Y3 = XofNs(NOCs(n, 3), 2);
      X4 = XofNs(NOCs(n, 4), 1);
      Y4 = XofNs(NOCs(n, 4), 2);

      SE=zeros(8);
      BB=zeros(3,8);

      for p = 1 :4
          % nodal coordinates
         
          %%%%%%%%%%%integration points
          XI = XNI(1,p);
          ETA = XNI(2,p);
          TJ11 = ((1 - ETA) * (X2 - X1) + (1 + ETA) * (X3 - X4)) / 4;
          TJ12 = ((1 - ETA) * (Y2 - Y1) + (1 + ETA) * (Y3 - Y4)) / 4;
          TJ21 = ((1 - XI) * (X4 - X1) + (1 + XI) * (X3 - X2)) / 4;
          TJ22 = ((1 - XI) * (Y4 - Y1) + (1 + XI) * (Y3 - Y2)) / 4;
          DJ = TJ11 * TJ22 - TJ12 * TJ21;
          %  --- Inversion of JACOBIAN
          A(1, 1) = TJ22 / DJ;
          A(1, 2) = -TJ12 / DJ;
          A(2, 1) = -TJ21 / DJ;
          A(2, 2) = TJ11 / DJ;
          %  --- Local derivatives of shape function
          g(1, 1) = -(1 - ETA) / 4;
          g(2, 1) = -(1 - XI) / 4;
          g(1, 2) = (1 - ETA) / 4;
          g(2, 2) = -(1 + XI) / 4;
          g(1, 3) = (1 + ETA) / 4;
          g(2, 3) = (1 + XI) / 4;
          g(1, 4) = -(1 + ETA) / 4;
          g(2, 4) = (1 - XI) / 4;
          
          B = A * g;
          %  --- BB = pose the B in the strainBB
          k=1:4;
          BB(1,2*k-1)=B(1,k);
          BB(3,2*k-1)=B(2,k);
          BB(2,2*k)=B(2,k);
          BB(3,2*k)=B(1,k);

          SE=SE+BB.'*D*BB*abs(DJ);
          %  --- strain_BB'*D= B'D
          end
end
     
           

