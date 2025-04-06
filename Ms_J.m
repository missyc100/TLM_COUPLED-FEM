function ME=Ms_J(n,NOCs,XofNs,Material_Elasstic,materils)
num_mater= materils(n);
den2=cell2mat(Material_Elasstic(num_mater,1));

XNI(1, 1) = -0.57735026919;
XNI(2, 1) = -0.57735026919;
XNI(1, 2) =  0.57735026919;
XNI(2, 2) = -0.57735026919;
XNI(1, 3) =  0.57735026919;
XNI(2, 3) =  0.57735026919;
XNI(1, 4) = -0.57735026919;
XNI(2, 4) =  0.57735026919;
% nodal coordinates
X1 = XofNs(NOCs(n, 1), 1);
Y1 = XofNs(NOCs(n, 1), 2);
X2 = XofNs(NOCs(n, 2), 1);
Y2 = XofNs(NOCs(n, 2), 2);
X3 = XofNs(NOCs(n, 3), 1);
Y3 = XofNs(NOCs(n, 3), 2);
X4 = XofNs(NOCs(n, 4), 1);
Y4 = XofNs(NOCs(n, 4), 2);

ME=zeros(8);
%integration points
for p = 1 :4
    XI = XNI(1,p);
    ETA = XNI(2,p);
    TJ11 = ((1 - ETA) * (X2 - X1) + (1 + ETA) * (X3 - X4)) / 4;
    TJ12 = ((1 - ETA) * (Y2 - Y1) + (1 + ETA) * (Y3 - Y4)) / 4;
    TJ21 = ((1 - XI) * (X4 - X1) + (1 + XI) * (X3 - X2)) / 4;
    TJ22 = ((1 - XI) * (Y4 - Y1) + (1 + XI) * (Y3 - Y2)) / 4;
    DJ = TJ11 * TJ22 - TJ12 * TJ21;
    ShapeFunM(1,1)=0.25*(1-XI)*(1-ETA);
    ShapeFunM(1,2)=0.25*(1+XI)*(1-ETA);
    ShapeFunM(1,3)=0.25*(1+XI)*(1+ETA);
    ShapeFunM(1,4)=0.25*(1-XI)*(1+ETA);
    NN1=[ShapeFunM(1,1)*eye(2), ShapeFunM(1,2)*eye(2), ShapeFunM(1,3)*eye(2), ShapeFunM(1,4)*eye(2)];
    ME=ME+den2*NN1.'*NN1*abs(DJ);
end
ME=diag(sum(ME,2));
end



