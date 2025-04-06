function [Me] = F_Me(n,Yy,materils_boundary,Material_Elasstic)
Y=Yy(:,2);
Me=zeros (4,4);
Y1= Y(n);
Y2= Y(n+1);
dy=(Y2-Y1);
num_mater=materils_boundary(n);
density=cell2mat(Material_Elasstic(num_mater,1));

Me(1,1)=2;    Me(1,2)=0;       Me(1,3)=1;      Me(1,4)=0;  
Me(2,1)=0;    Me(2,2)=2;       Me(2,3)=0;      Me(2,4)=1;  
Me(3,1)=1;    Me(3,2)=0;       Me(3,3)=2;      Me(3,4)=0;  
Me(4,1)=0;    Me(4,2)=1;       Me(4,3)=0;      Me(4,4)=2;                 
Me=Me*density*dy/6;   
end

