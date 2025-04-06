%%              Coupled TLM with halfspace 阻抗
clc
clear
% 内域几何信息
NNs=121;
NEs=100; % 从inp 文件查看结点号和单元号
NOCs=dlmread('element_Solid.txt',',',[1,1,NEs,4]);
XofNs=dlmread('node_Solid.txt',',',[1,1,NNs,2]);

% 材料信息
cp=[766;766];cs=[442;442];rs=[2700;2700];
Type_material=size(rs,1);
lamds=(-2*cs.^2+cp.^2).*rs;
Gs=cs.^2.*rs;
vt=lamds./(2*(lamds+Gs));
dims=2; % 计算维数
[Material_Elasstic] = Materials_matrix(rs,vt,cs,Type_material,dims);
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

%                       内域刚度矩阵计算
[MM1,KK1] = integrate_Solid(NOCs,XofNs,Material_Elasstic,materils);  %J denotes the lumped
CC1=sparse(size(MM1,1),size(MM1,1),0);

% 基岩强制位移
ntime=2001; dt=0.01;
t0=2;fp=0.7;A_max=0.001;
t=(0:dt:(ntime-1)*dt);
disp=A_max*(1-2*(pi*fp*(t-t0)).^2).*exp(-(pi*fp*(t-t0)).^2);
disp(abs(disp)<1.0e-8)=0;
incidence=fft(disp);
Tp=(ntime-1)*dt;
d_omiga=2*pi/Tp;
omiga=0:d_omiga:(ntime-1)*d_omiga;
omiga(1)=1.0e-10;
Feq=zeros(2*NNs,ntime);
plot(omiga,abs(incidence))
fix_finitedomaiin= importdata('Bedrock_node.txt');
fix_finitedomaiin= fix_finitedomaiin(find( ~ isnan(fix_finitedomaiin)));
fix_finitedomaiin=[fix_finitedomaiin(1):fix_finitedomaiin(end):fix_finitedomaiin(2)].';
fix_fitedomaiX=fix_finitedomaiin*2-1;
fix_fitedomaiY=fix_finitedomaiin*2;

UU=zeros(size(MM1,1),ntime);
%% TLM
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
[AAL1,ABL1,AGL,AML,ADL1] = Integrate_TLM(BNL,NOC_BL,Xof_B1,materils_boundary,Material_Elasstic);
ABL1=1i*ABL1;
% 刚度矩阵
SBL=eye(2*BNL,2*BNL);
DofL_exclude_fix=2*BNL-2;
fix_xL=BNL*2-1;
fix_yL=BNL*2;
fix_L=[fix_xL;fix_yL];
% 自由场
Free_FL=zeros(2*BNL,ntime);
Free_UL=zeros(2*BNL,ntime);
Free_ForceL=zeros(2*BNL,ntime);

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
[AAR1,ABR1,AGR,AMR,ADR1] = Integrate_TLM(BNR,NOC_BR,Xof_B2,materils_boundary,Material_Elasstic);
ABR1=1i*ABR1;
% 刚度矩阵
SBR=eye(2*BNR,2*BNR);
DofR_exclude_fix=2*BNR-2;
esp2=1.0e-8;
fix_xR=BNR*2-1;
fix_yR=BNR*2;
fix_R=[fix_xR;fix_yR];
% 自由场
Free_FR=zeros(2*BNR,ntime);
Free_UR=zeros(2*BNR,ntime);
Free_ForceR=zeros(2*BNR,ntime);


AAL1(fix_L,:)=[]; AAL1(:,fix_L)=[];
ABL1(fix_L,:)=[]; ABL1(:,fix_L)=[];
AAR1(fix_R,:)=[]; AAR1(:,fix_R)=[];
ABR1(fix_R,:)=[]; ABR1(:,fix_R)=[];

for jj=1:round(ntime/2)+1
    jj
    Fre=omiga(jj);
    %%———————————————————————————左右侧自由场求解
    % 左侧自由场计算
    ACL=AGL-Fre^2*AML;
    Free_FL(:,jj)=Free_FL(:,jj)-ACL(:,fix_xL)*incidence(jj);
    ACL(fix_xL,:)=0;  %
    ACL(:,fix_xL)=0;
    ACL(fix_xL,fix_xL)=1;
    Free_FL(fix_xL,jj)=incidence(jj);
    ACL(fix_yL,:)=0;
    ACL(:,fix_yL)=0;
    ACL(fix_yL,fix_yL)=1;
    Free_FL(fix_yL,jj)=0;
    Free_UL(:,jj)=ACL\Free_FL(:,jj);
    Free_ForceL(:,jj)=ADL1*Free_UL(:,jj); 
    ADL =ADL1;
    % 右侧自由场计算
    ACR=AGR-Fre^2*AMR;
    Free_FR(:,jj)=Free_FR(:,jj)-ACR(:,fix_xR)*incidence(jj);
    ACR(fix_xR,:)=0;
    ACR(:,fix_xR)=0;
    ACR(fix_xR,fix_xR)=1;
    Free_FR(fix_xR,jj)=incidence(jj);
    ACR(fix_yR,:)=0;
    ACR(:,fix_yR)=0;
    ACR(fix_yR,fix_yR)=1;
    Free_FR(fix_yR,jj)=0;
    Free_UR(:,jj)=ACR\Free_FR(:,jj);
    Free_ForceR(:,jj)=ADR1*Free_UR(:,jj); 
    %%———————————————————————— 无限域动力刚度矩阵计算
    %% 左侧

    ACL(fix_L,:)=[]; ACL(:,fix_L)=[];
    ADL(fix_L,:)=[]; ADL(:,fix_L)=[];
    %%%%%%复模态%%%%%%%%%%
    MB(1:DofL_exclude_fix,1:DofL_exclude_fix)=-ACL;%%%%%
    MB(DofL_exclude_fix+1:2*DofL_exclude_fix,DofL_exclude_fix+1:2*DofL_exclude_fix)=AAL1;
    
    KB(1:DofL_exclude_fix,DofL_exclude_fix+1:2*DofL_exclude_fix)=-ACL;
    KB(DofL_exclude_fix+1:2*DofL_exclude_fix,1:DofL_exclude_fix)=-ACL;
    KB(DofL_exclude_fix+1:2*DofL_exclude_fix,DofL_exclude_fix+1:2*DofL_exclude_fix)=-ABL1;
    Ae=MB\KB;
    [Fei,lad]=eig(Ae); %每一频率下的特征向量矩阵、和特征值矩阵
    % 求左边动力刚度矩阵
    jj4=0;
    for jj2=1:2*2*(BNL-1)%%%%%%%找出模态中虚特征值小于0对应的的特征向量
        if  ((abs(imag(lad(jj2,jj2)))<esp2) && (real((lad(jj2,jj2)))<-esp2 ))...
                || (imag(lad(jj2,jj2))>esp2)
            jj4=1+jj4;
            LaL(jj4,jj4)=lad(jj2,jj2);
            Fb_L(:,jj4)=Fei(:,jj2);%%%%%%%%%%找出特征值虚部小于零的特征向量
        end
    end
    Number_jj1(jj)=jj4;
    FbL=Fb_L(1:2*(BNL-1),:);%
    Sb_L=-AAL1*FbL*(1i)*LaL/(FbL)-ADL; %
    SBL(1:end-2,1:end-2)=Sb_L;

    %% 右侧
    ADR =ADR1;
    ACR(fix_R,:)=[]; ACR(:,fix_R)=[];
    ADR(fix_R,:)=[]; ADR(:,fix_R)=[];
    %%%%%%复模态%%%%%%%%%%
    MBR(1:DofR_exclude_fix,1:DofR_exclude_fix)=-ACR;%%%%%
    MBR(DofR_exclude_fix+1:2*DofR_exclude_fix,DofR_exclude_fix+1:2*DofR_exclude_fix)=AAR1;
    
    KBR(1:DofR_exclude_fix,DofR_exclude_fix+1:2*DofR_exclude_fix)=-ACR;
    KBR(DofR_exclude_fix+1:2*DofR_exclude_fix,1:DofR_exclude_fix)=-ACR;
    KBR(DofR_exclude_fix+1:2*DofR_exclude_fix,DofR_exclude_fix+1:2*DofR_exclude_fix)=-ABR1;
    AeR=MBR\KBR;
    %%%%%%%%%%%求每一频率下的特征值和特征向量
    [FeiR,ladR]=eig(AeR);
    %%%%%%%%%%%%%每一频率下的特征向量矩阵、和特征值矩阵
    % 求左边动力刚度矩阵
    jj4=0;
    for jj2=1:2*2*(BNR-1)%%%%%%%找出模态中虚特征值小于0对应的的特征向量
        if  ((abs(imag(ladR(jj2,jj2)))<esp2) && (real((ladR(jj2,jj2)))>esp2 ))...
                || (imag(ladR(jj2,jj2))<-esp2)
            jj4=1+jj4;
            LaR(jj4,jj4)=ladR(jj2,jj2);
            Fb_R(:,jj4)=FeiR(:,jj2);%%%%%%%%%%找出特征值虚部小于零的特征向量
        end
    end
    Number_jj1(jj)=jj4;
    FbR=Fb_R(1:2*(BNR-1),:);%
    Sb_R=AAR1*FbR*(1i)*LaR/(FbR)+ADR; %
    SBR(1:end-2,1:end-2)=Sb_R;

  %% 有限域+边界
    [Keq]=Assembly_Finite_infinite(Fre,MM1,CC1,KK1,SBL,SBR,plus_BLL,plus_BRR);
    Feq(plus_BRR,jj)=SBR*Free_UR(:,jj)-Free_ForceR(:,jj);
    Feq(plus_BLL,jj)=SBL*Free_UL(:,jj)+Free_ForceL(:,jj);

     for i=1:length(fix_fitedomaiX)
         Feq(:,jj)=Feq(:,jj)-Keq(:,fix_fitedomaiX(i))*incidence(jj);
         Keq(fix_fitedomaiX(i),:)=0;
         Keq(:,fix_fitedomaiX(i))=0;
         Keq(fix_fitedomaiX(i),fix_fitedomaiX(i))=1;
         Feq(fix_fitedomaiX(i),jj)=incidence(jj);

         Keq(fix_fitedomaiY(i),:)=0;
         Keq(:,fix_fitedomaiY(i))=0;
         Keq(fix_fitedomaiY(i),fix_fitedomaiY(i))=1;
         Feq(fix_fitedomaiY(i),jj)=0;
    end
    [U1,isornot] = solve(sparse(Keq),sparse(Feq(:,jj)));
    UU(:,jj)=U1;
end

reference=[1:11].';
result=UU(2*reference-1,:).';
result=ifft(result,'symmetric');

% 自由场
freefield_result=ifft(Free_UR.','symmetric');

plot(t.',freefield_result(:,1))
hold on
plot(t.',result(:,:))
hold off
legend('自由场','有限元')


