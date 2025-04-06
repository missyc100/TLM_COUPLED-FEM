function [D2] = F_D2(n )

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
D2=zeros (4,4);                                                                                             
D2(1,1)=G;        D2(1,2)=0;         D2(1,3)=-G;             D2(1,4)=0;         
D2(2,1)=0;        D2(2,2)=A+2*G;     D2(2,3)=0;              D2(2,4)=-(A+2*G);  
D2(3,1)=G;        D2(3,2)=0;         D2(3,3)=-G;             D2(3,4)=0;        
D2(4,1)=0;        D2(4,2)=A+2*G;     D2(4,3)=0;              D2(4,4)=-(A+2*G);                  
D2=-ele_nx/ele_ny*D2/2;   
end

