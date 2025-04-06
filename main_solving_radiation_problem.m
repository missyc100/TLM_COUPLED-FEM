%%                               TLM coupled FEM
clc
clear
%% 内域几何信息
NNs=121;
NEs=100; % 从inp 文件查看结点号和单元号
NOCs=dlmread('element_Solid.txt',',',[1,1,NEs,4]);
XofNs=dlmread('node_Solid.txt',',',[1,1,NNs,2]);

%% 材料信息
cp=[400;400];cs=[250;250];rs=[2700;2700];
Type_material=size(rs,1);
lamds=(-2*cs.^2+cp.^2).*rs;
Gs=cs.^2.*rs;
vt=lamds./(2*(lamds+Gs));
[Material_Elasstic] = Materials_matrix(rs,vt,cs,Type_material);
num_mater=[1:Type_material].';
materils=zeros(NEs,1);
  % 各个单元材料属性
ele_soil1= importdata('element_soil1.txt');
ele_soil1= ele_soil1(find( ~ isnan(ele_soil1)));
ele_soil1=[ele_soil1(1):ele_soil1(end):ele_soil1(2)].';
materils(ele_soil1,1)=num_mater(1);
  % 土层类单元
ele_soil2 = setdiff([1:NEs], [ele_soil1;]);
materils(ele_soil2,1)=num_mater(2);

%%                        内域刚度矩阵计算
[MM1,KK1] = integrate_Solid(NOCs,XofNs,Material_Elasstic,materils);  %J denotes the lumped
 CC1=sparse(size(MM1,1),size(MM1,1),0);


%%                           所加的外荷载波源信息
% 1. dirac 脉冲
% T =5;
% fp=1/T;
% x = t./T;
% Force = A_max*16.*(x.^3.*heaviside(x)-4.*(x-0.25).^3.*heaviside(x-0.25)+6.*(x-0.5).^3.*heaviside(x-0.5)-4.*(x-0.75).^3.*heaviside(x-0.75)+(x-1).^3.*heaviside(x-1));
% Force((abs(Force)<1.0e-8))=0;
% 2. Ricker 小波
ntime=1024*32; dt=0.01;
t0=4;fp=0.5;A_max=100000;
t=(0:dt:(ntime-1)*dt);
Force=A_max*(1-2*(pi*fp*(t-t0)).^2).*exp(-(pi*fp*(t-t0)).^2);
Force(abs(Force)<1.0e-8)=0;

% 
incidence=fft(Force.');
Tp=(ntime-1)*dt;
d_omiga=2*pi/Tp;
omiga=0:d_omiga:(ntime-1)*d_omiga;
omiga(1)=1.0e-8;
plot(omiga./(2*pi),abs(incidence));
% 加载点
point_force=[100,0];
Dis=pdist2(XofNs, point_force,'euclidean');
[minDist, node_num_force] = min(Dis);
UU=zeros(size(MM1,1),ntime);
Feq=zeros(size(MM1,1),ntime);
Feq(2*node_num_force,:)=ones(size(node_num_force,1),1)*incidence;
% 3. 底部基岩固定边界条件
fix_finitedomain= importdata('Bedrock_node.txt');
fix_finitedomain= fix_finitedomain(find( ~ isnan(fix_finitedomain)));
fix_finitedomain=[fix_finitedomain(1):fix_finitedomain(end):fix_finitedomain(2)].';
fix_fitedomain=[fix_finitedomain*2-1;fix_finitedomain*2];

%%                             TLM
%% ——————左侧信息
N_B1= importdata('BLL_node.txt');
N_B1= N_B1(find( ~ isnan(N_B1)));
N_B1=[N_B1(1):N_B1(end):N_B1(2)].';
BNL=size(N_B1,1);
Xof_B1=XofNs(N_B1,:);
NOC_BL=[1:BNL-1;2:BNL]';
NEL=size(NOC_BL,1);
B1_inf=[N_B1,Xof_B1];
a_B1_inf=sortrows(B1_inf,3);
plus_BLL=zeros(2*BNL,1);
plus_BLL(1:2:end,1)=2*a_B1_inf(:,1)-1;
plus_BLL(2:2:end,:)=2*a_B1_inf(:,1);
%_边界外围土层深度；
materils_boundary=zeros(NEL,1);
soil_layer1=[60];  %60为土层1深度
[a,node_soil_layer1]=ismember(soil_layer1,a_B1_inf(:,3));
ele_soil1=[1:node_soil_layer1(1)-1].';
ele_soil2=setdiff([1:NEL].',ele_soil1);
materils_boundary(ele_soil1,1)=num_mater(1);
materils_boundary([ele_soil2],1)=num_mater(2);
%____TLM左侧动力刚度相关矩阵
[AAL1,ABL1,AGL1,AML1,ADL1] = Integrate_TLM(BNL,NOC_BL,Xof_B1,materils_boundary,Material_Elasstic);
ABL1=1i*ABL1;
% 刚度矩阵
SBL=eye(2*BNL,2*BNL);
DofL_exclude_fix=2*BNL-2;
fix_xL=BNL*2-1;
fix_yL=BNL*2;
fix_L=[fix_xL;fix_yL];
Z1=zeros(DofL_exclude_fix,DofL_exclude_fix);
%% —————TLM右侧动力刚度相关矩阵
N_B2= importdata('BRR_node.txt');
N_B2= N_B2(find( ~ isnan(N_B2)));
N_B2=[N_B2(1):N_B2(end):N_B2(2)].';
BNR=size(N_B2,1);
Xof_B2=XofNs(N_B2,:);
NOC_BR=[1:BNR-1;2:BNR]';
NER=size(NOC_BR,1);
B2_inf=[N_B2,Xof_B2];
a_B2_inf=sortrows(B2_inf,3);
plus_BRR=zeros(2*BNR,1);
plus_BRR(1:2:end,1)=2*a_B2_inf(:,1)-1;
plus_BRR(2:2:end,:)=2*a_B2_inf(:,1);

materils_boundary=zeros(NER,1);
soil_layer1=[60];  %60为土层1深度
[a,node_soil_layer1]=ismember(soil_layer1,a_B2_inf(:,3));
ele_soil1=[1:node_soil_layer1(1)-1].';
ele_soil2=setdiff([1:NER].',ele_soil1);
materils_boundary(ele_soil1,1)=num_mater(1);
materils_boundary([ele_soil2],1)=num_mater(2);
%____TLM右侧侧动力刚度相关矩阵
[AAR1,ABR1,AGR1,AMR1,ADR1] = Integrate_TLM(BNR,NOC_BR,Xof_B2,materils_boundary,Material_Elasstic);
ABR1=1i*ABR1;
% 刚度矩阵4
SBR=eye(2*BNR,2*BNR);
DofR_exclude_fix=2*BNR-2;
esp2=1.0e-8;
fix_xR=BNR*2-1;
fix_yR=BNR*2;
fix_R=[fix_xR;fix_yR];
Zr=zeros(DofR_exclude_fix,DofR_exclude_fix);

%    剔除左右侧边界固定节点号
AAL1(fix_L,:)=[]; AAL1(:,fix_L)=[];
ABL1(fix_L,:)=[]; ABL1(:,fix_L)=[];
AGL1(fix_L,:)=[]; AGL1(:,fix_L)=[];
AML1(fix_L,:)=[]; AML1(:,fix_L)=[];
ADL1(fix_L,:)=[]; ADL1(:,fix_L)=[];

AAR1(fix_R,:)=[]; AAR1(:,fix_R)=[];
ABR1(fix_R,:)=[]; ABR1(:,fix_R)=[];
AGR1(fix_R,:)=[]; AGR1(:,fix_R)=[];
AMR1(fix_R,:)=[]; AMR1(:,fix_R)=[];
ADR1(fix_R,:)=[]; ADR1(:,fix_R)=[];
%%            
for jj=1:ntime/2+1
    jj
    Fre=omiga(jj);
    ACL1=AGL1-Fre^2*AML1;
    %%———————————————————————— 无限域动力刚度矩阵计算
    %% 左侧复模态%%%%%%%%%%
    MB=[-ACL1,Z1;Z1,AAL1];
    KB=[Z1,-ACL1;-ACL1,-ABL1];
    Ae=MB\KB;
    [Fei,lad]=eig(Ae); %每一频率下的特征向量矩阵、和特征值矩阵
    % 求左边动力刚度矩阵
    [LaL,FbL,jj4]=solve_eig_L(Fei,lad,DofL_exclude_fix,esp2);
    Sb_L=AAL1*FbL*(1i)*LaL/(FbL)+ADL1;
    SBL(1:end-2,1:end-2)=Sb_L;

    %% 右侧复模态%%%%%%%%%%
    ACR1=AGR1-Fre^2*AMR1;
    MBR=[-ACR1,Zr;Zr,AAR1];
    KBR=[Zr,-ACR1;-ACR1,-ABR1];
    AeR=MBR\KBR;
    %求每一频率下的特征值和特征向量
    [FeiR,ladR]=eig(AeR);
    [LaR,FbR,jj4R]=solve_eig_R(FeiR,ladR,DofR_exclude_fix,esp2);
    Sb_R=AAR1*FbR*(1i)*LaR/(FbR)+ADR1; %
    SBR(1:end-2,1:end-2)=Sb_R;

    %% 有限域+边界
    [Keq]=Assembly_Finite_infinite(Fre,MM1,CC1,KK1,-SBL,SBR,plus_BLL,plus_BRR);

    Keq(fix_fitedomain,:)=0;
    Keq(:,fix_fitedomain)=0;
    Keq(fix_fitedomain,fix_fitedomain)=eye(size(fix_fitedomain,1));
    Feq(fix_fitedomain,jj)=0;

    % [U1,isornot] = solve(sparse(Keq),sparse(Feq(:,jj)));
    % UU(:,jj)=U1;
    UU(:,jj)=Keq\Feq(:,jj);
    % UU(:,jj)=U1;
end

reference=[1:11].';
result=UU(2*reference,:).';
result=ifft(result,'symmetric');




