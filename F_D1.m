function [D1] = F_D1(n,Yy,materils_boundary,Material_Elasstic) 
num_mater=materils_boundary(n);
A=cell2mat(Material_Elasstic(num_mater,2));
G=cell2mat(Material_Elasstic(num_mater,3));

D1=zeros (4,4);                                                                                             
D1(1,1)=0;        D1(1,2)=-A;         D1(1,3)=0;             D1(1,4)=A;         
D1(2,1)=-G;       D1(2,2)=0;          D1(2,3)=G;             D1(2,4)=0;  
D1(3,1)=0;        D1(3,2)=-A;         D1(3,3)=0;             D1(3,4)=A;        
D1(4,1)=-G;       D1(4,2)=0;          D1(4,3)=G;             D1(4,4)=0;                  
D1=-D1/2;   
end

