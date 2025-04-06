function [A1] =F_A1(n,Yy,materils_boundary,Material_Elasstic)
Y=Yy(:,2);
A1=zeros (4,4);
Y1= Y(n);Y2= Y(n+1);
dy=(Y2-Y1);
num_mater=materils_boundary(n);
A=cell2mat(Material_Elasstic(num_mater,2));
G=cell2mat(Material_Elasstic(num_mater,3));
A1(1,1)=2*(A+2*G);       A1(1,2)=0;         A1(1,3)=A+2*G;    A1(1,4)=0;
A1(2,1)=0;               A1(2,2)=2*G;       A1(2,3)=0;        A1(2,4)=G;
A1(3,1)=A+2*G;           A1(3,2)=0;         A1(3,3)=2*(A+2*G);A1(3,4)=0;
A1(4,1)=0;               A1(4,2)=G;         A1(4,3)=0;        A1(4,4)=2*G;
A1=A1*dy/6;
end

