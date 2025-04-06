function [A3] =F_A3(n)
global   nxL a_B1_inf  solid_layer1 solid_layer2
Y=a_B1_inf(:,3);
A3=zeros (4,4);
Y1= Y(n);Y2= Y(n+1);
dy=(Y2-Y1);
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
A3(1,1)=0;         A3(1,2)=1;       A3(1,3)=0;     A3(1,4)=1;
A3(2,1)=1;         A3(2,2)=0;       A3(2,3)=1;     A3(2,4)=0;
A3(3,1)=0;         A3(3,2)=1;       A3(3,3)=0;     A3(3,4)=1;
A3(4,1)=1;         A3(4,2)=0;       A3(4,3)=1;     A3(4,4)=0;
A3=-(A+G)*(ele_nx/ele_ny)*A3*dy/4;
end

