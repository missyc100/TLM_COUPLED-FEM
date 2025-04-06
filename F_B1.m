function [B1] = F_B1(n,Yy,materils_boundary,Material_Elasstic)
num_mater=materils_boundary(n);
A=cell2mat(Material_Elasstic(num_mater,2));
G=cell2mat(Material_Elasstic(num_mater,3));
B1=zeros(4,4);
B1(1,1)=0;         B1(1,2)=(A-G);     B1(1,3)=0;        B1(1,4)=-(A+G);  
B1(2,1)=-(A-G);    B1(2,2)=0;       B1(2,3)=-(A+G);   B1(2,4)=0;        
B1(3,1)=0;         B1(3,2)=(A+G);   B1(3,3)=0;        B1(3,4)=-(A-G);  
B1(4,1)=(A+G);     B1(4,2)=0;       B1(4,3)=(A-G);    B1(4,4)=0;        
B1=-B1/2;   
end

