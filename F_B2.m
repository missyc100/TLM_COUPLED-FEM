function [B2] = F_B2(n )
global  nxL  a_B1_inf   solid_layer1 solid_layer2
Y=a_B1_inf(:,3);
Y1= Y(n);Y2= Y(n+1);
mid_Y=(Y1+Y2)/2;
if mid_Y<=solid_layer1
    A=dlmread('Solid_materials.txt',',',[1,0,1,0]);
    G=dlmread('Solid_materials.txt',',',[1,1,1,1]);
elseif (mid_Y>solid_layer1)  && (mid_Y<solid_layer2)
    A=dlmread('Solid_materials.txt',',',[2,0,2,0]);
    G=dlmread('Solid_materials.txt',',',[2,1,2,1]);
else
    A=dlmread('Solid_materials.txt',',',[3,0,3,0]);
    G=dlmread('Solid_materials.txt',',',[3,1,3,1]);
end
ele_nx=nxL(n,1); ele_ny=nxL(n,2);
B2=zeros(4,4);
B2(1,1)=0;         B2(1,2)=0;         B2(1,3)=G;    B2(1,4)=0;  
B2(2,1)=0;         B2(2,2)=0;         B2(2,3)=0;    B2(2,4)=(A+2*G);        
B2(3,1)=-G;        B2(3,2)=0;         B2(3,3)=0;    B2(3,4)=0;  
B2(4,1)=0;         B2(4,2)=-(A+2*G);  B2(4,3)=0;    B2(4,4)=0;        
B2=-(ele_nx/ele_ny)*B2;   
end

