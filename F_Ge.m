function [Ge] = F_Ge(n,Yy,materils_boundary,Material_Elasstic)

Y=Yy(:,2);
Ge=zeros (4,4);
Y1= Y(n);Y2= Y(n+1);
dy=(Y2-Y1);
num_mater=materils_boundary(n);
G=cell2mat(Material_Elasstic(num_mater,3));
A=cell2mat(Material_Elasstic(num_mater,2));

Ge(1,1)=G;        Ge(1,2)=0;         Ge(1,3)=-G;    Ge(1,4)=0;            
Ge(2,1)=0;        Ge(2,2)=A+2*G;     Ge(2,3)=0;     Ge(2,4)=-(A+2*G);  
Ge(3,1)=-G;       Ge(3,2)=0;         Ge(3,3)=G;     Ge(3,4)=0; 
Ge(4,1)=0;        Ge(4,2)=-(A+2*G);  Ge(4,3)=0;     Ge(4,4)=A+2*G;                  
Ge=Ge/dy;   
end

