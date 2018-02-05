function [ Motifout ] = Position_motifs(M,V_motifs)
% Count the frequency with which each node occurs in positions within 
% bipartite motifs up to six nodes
%   M is a binary two-dimensional matrix
%   V_motifs: is a vector with the id of motif positions you want count (between 1
%   and 148)
% Output: a table with one row for each motif position and a number of columns 
% equal to 1 + the sum of the dimensions of M (number of rows + number of columns). 
% The first column (IDmotif) gives the ID of the motif position, 
% the Ax columns give the frequency with which each row node occurs in each
% motif position, the Bx columns give the frequency with which each column 
% node occurs in each motif position. NaN is returned when a node cannot 
% occur in a given position because that position can only be occupied by 
% nodes in the other level. For example, position 1 (IDmotif 1) can only be
% occupied by column nodes and therefore a count of NaN is returned for all 
% row nodes for this position

if nargin==1
    if Check_pos(M)==0
        return
    end
elseif nargin==2
    if Check_pos(M,V_motifs)==0
        return
    end
end

Mot_var_id={};
Mot_var_id{1}='IDmotif';

for k=1:size(M,1)
    Mot_var_id{k+1}=strcat('A',num2str(k));
end
lab_p={};
for k=1:size(M,2)
    Mot_var_id{k+size(M,1)+1}=strcat('B',num2str(k));
end



mat=struct;
mat.M=M;

mat.z=size(M,1);
mat.znul=nan(size(M,1),1);
mat.p=size(M,2);
mat.pnul=nan(size(M,2),1);
mat.jP=ones(mat.p,1);
mat.jZ=ones(mat.z,1);
mat.J=ones(mat.z,mat.p);
mat.JP=ones(mat.p,mat.p);
mat.JZ=ones(mat.z,mat.z);

mat.MT=mat.M';
mat.N=mat.J-mat.M;
mat.NT=mat.N';


mat.dZ=mat.M*mat.jP;
mat.dP=mat.MT*mat.jZ;

mat.Z=mat.M*mat.MT;
mat.Y=mat.M*mat.NT;
mat.X=mat.N*mat.MT;

mat.P=mat.MT*mat.M;
mat.Q=mat.MT*mat.N;
mat.R=mat.NT*mat.M;

%%
mat.JP3=ones(mat.p,mat.p,mat.p);
mat.JZ3=ones(mat.z,mat.z,mat.z);
mat.KP3=mat.JP3;
for (i=1 : mat.p)
    for (j=1 : mat.p)
        mat.KP3(i,j,j)=0;
        mat.KP3(j,i,j)=0;
        mat.KP3(j,j,i)= 0;
    end
end
mat.KZ3=mat.JZ3
for (i =1 : mat.z)
    for (j=1 : mat.z)
        mat.KZ3(i,j,j)=0;
        mat.KZ3(j,i,j)=0;
        mat.KZ3(j,j,i)=0;
    end
end



mat.AZ=tensor_make(mat.M,mat.M);
mat.BZ=tensor_make(mat.M,mat.N);
mat.CZ=tensor_make(mat.N,mat.M);
mat.DZ=tensor_make(mat.N,mat.N);

mat.AP=tensor_make(mat.MT,mat.MT);
mat.BP=tensor_make(mat.MT,mat.NT);
mat.CP=tensor_make(mat.NT,mat.MT);
mat.DP=tensor_make(mat.NT,mat.NT);


mat.MTA = tensorR(mat.MT,mat.AZ);
mat.MTB = tensorR(mat.MT,mat.BZ);
mat.MTC = tensorR(mat.MT,mat.CZ);
mat.MTD = tensorR(mat.MT,mat.DZ);

mat.MA = tensorR(mat.M,mat.AP);
mat.MB = tensorR(mat.M,mat.BP);
mat.MC = tensorR(mat.M,mat.CP);
mat.MD = tensorR(mat.M,mat.DP);

mat.NTA = tensorR(mat.NT,mat.AZ);
mat.NTB = tensorR(mat.NT,mat.BZ);
mat.NTC = tensorR(mat.NT,mat.CZ);

mat.Na = tensorR(mat.N,mat.AP);
mat.NB = tensorR(mat.N,mat.BP);
mat.NC = tensorR(mat.N,mat.CP);



count=1;
for n_motifs=V_motifs
    Motif_r=0;
    Motif_int=0;
    Motif=0;
    %%Motifs 1 count
    
    fh=str2func(strcat('Motifs_pos',num2str(n_motifs)));
    Motif=fh(mat);
    
    
    Motifout(count,:)=[n_motifs,Motif];
    count=count+1;
end

Motifout=array2table(Motifout,'VariableNames',Mot_var_id);

%Motifout=array2table(Motifv)
end


function [ Mot ] = Motifs_pos1(mat)
Mot=[ mat.znul' mat.dP'];
end

function [ Mot ] = Motifs_pos2(mat)
Mot=[mat.dZ' mat.pnul'];
end

function [ Mot ] = Motifs_pos3(mat)
Mot= mat.P*mat.jP - mat.dP;
Mot=[mat.znul' Mot'];
end

function [ Mot ] = Motifs_pos4(mat)
Mot= mat.dZ.* (mat.dZ - mat.jZ) / 2;
Mot=[ Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos5(mat)
Mot= mat.dP.* (mat.dP - mat.jP) / 2;
Mot=[ mat.znul' Mot' ];
end

function [ Mot ] = Motifs_pos6(mat)
Mot= mat.Z*mat.jZ - mat.dZ;
Mot=[  Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos7(mat)
Mot= mat.dP .* (mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP) / 6;
Mot=[  mat.znul' Mot' ];
end

function [ Mot ] = Motifs_pos8(mat)
Mot= mat.M * ((mat.dP - mat.jP) .* (mat.dP - 2 * mat.jP)) / 2;
Mot=[  Mot'  mat.pnul'  ];
end

function [ Mot ] = Motifs_pos9(mat)
Mot=(mat.P .* mat.R) * mat.jP;
Mot=[  mat.znul' Mot'    ];
end

function [ Mot ] = Motifs_pos10(mat)
Mot=(mat.P .* mat.Q) * mat.jP;
Mot=[ mat.znul'  Mot'    ];
end

function [ Mot ] = Motifs_pos11(mat)
Mot=(mat.X .* mat.Z) * mat.jZ;
Mot=[ Mot'    mat.pnul'  ];
end

function [ Mot ] = Motifs_pos12(mat)
Mot=(mat.Y .* mat.Z) * mat.jZ;
Mot=[ Mot'    mat.pnul'  ];
end

function [ Mot ] = Motifs_pos13(mat)
Mot=(mat.P .* (mat.P - mat.JP)) * mat.jP / 2 - mat.dP .* (mat.dP - mat.jP) / 2;
Mot=[ mat.znul' Mot'     ];
end

function [ Mot ] = Motifs_pos14(mat)
Mot=(mat.Z .* (mat.Z - mat.JZ)) * mat.jZ / 2 - mat.dZ .* (mat.dZ - mat.jZ) / 2;
Mot=[ Mot'    mat.pnul'  ];
end

function [ Mot ] = Motifs_pos15(mat)
Mot=mat.MT * ((mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ)) / 2;
Mot=[ mat.znul' Mot'     ];
end

function [ Mot ] = Motifs_pos16(mat)
Mot=mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) / 6;
Mot=[ Mot'    mat.pnul'  ];
end

function [ Mot ] = Motifs_pos17(mat)
Mot=mat.dP .* (mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP) .* (mat.dP - 3 .* mat.jP) / 24;
Mot=[ mat.znul'  Mot'     ];
end

function [ Mot ] = Motifs_pos18(mat)
Mot=mat.M * ((mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP) .* (mat.dP - 3 * mat.jP)) / 6;
Mot=[  Mot'    mat.pnul'  ];
end

function [ Mot ] = Motifs_pos19(mat)
Mot=(mat.P .* mat.R .* (mat.R - mat.JP)) * mat.jP / 2;
Mot=[  mat.znul'  Mot'    ];
end

function [ Mot ] = Motifs_pos20(mat)
Mot=(mat.P .* mat.Q .* (mat.Q - mat.JP)) * mat.jP / 2;
Mot=[  mat.znul'  Mot'    ];
end

function [ Mot ] = Motifs_pos21(mat)
Mot=mat.N .* (mat.M * ((mat.Q - mat.JP) .* mat.P)) * mat.jP;
Mot=[   Mot'   mat.pnul'  ];
end

function [ Mot ] = Motifs_pos22(mat)
Mot=(mat.M .* (mat.M * (mat.Q .* (mat.Q - mat.JP)))) * mat.jP / 2;
Mot=[   Mot'   mat.pnul'  ];
end

function [ Mot ] = Motifs_pos23(mat)
Mot=(mat.P .* mat.Q .* mat.R) * mat.jP;
Mot=[   mat.znul'  Mot'   ];
end

function [ Mot ] = Motifs_pos24(mat)
Mot=(mat.N .* (mat.M * (mat.P .* mat.R))) * mat.jP;
Mot=[Mot' mat.pnul'   ];
end

function [ Mot ] = Motifs_pos25(mat)
Mot=(mat.M .* (mat.M * (mat.Q .* mat.R))) * mat.jP / 2;
Mot=[Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos26(mat)
Mot=(mat.P .* (mat.P - mat.JP) .* mat.R) * mat.jP / 2;
Mot=[mat.znul' Mot'];
end

function [ Mot ] = Motifs_pos27(mat)
Mot=(mat.P .* (mat.P - mat.JP) .* mat.Q) * mat.jP / 2;
Mot=[mat.znul' Mot'];
end

function [ Mot ] = Motifs_pos28(mat)
Mot=(mat.N .* (mat.M * (mat.P .* (mat.P - mat.JP)))) * mat.jP / 2;
Mot=[ Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos29(mat)
Mot= (mat.M .* (mat.M * (mat.Q .* (mat.P - mat.JP)))) * mat.jP;
Mot=[ Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos30(mat)
Mot= (mat.P .* (mat.P - mat.JP) .* (mat.P - 2 * mat.JP)) * mat.jP / 6 -  mat.dP .* (mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP) / 6;
Mot=[ mat.znul' Mot' ];
end

function [ Mot ] = Motifs_pos31(mat)
Mot=(mat.M .* (mat.M * ((mat.P - mat.JP) .* (mat.P - 2 .* mat.JP))) * mat.jP - mat.M * ((mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP))) / 4;
Mot=[ Mot'  mat.pnul' ];
end


function [ Mot ] = Motifs_pos32(mat)
Mot=(mat.NT .* (mat.MT * ((mat.Y - mat.JZ) .* mat.Z))) * mat.jZ;
Mot=[mat.znul'  Mot'];
end


function [ Mot ] = Motifs_pos33(mat)
Mot=(mat.MT .* (mat.MT * (mat.Y .* (mat.Y - mat.JZ)))) * mat.jZ / 2;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos34(mat)
Mot=(mat.Z .* mat.X .* (mat.X - mat.JZ)) * mat.jZ / 2;
Mot=[ Mot' mat.pnul' ];
end

function [ Mot ] = Motifs_pos35(mat)
Mot=(mat.Z .* mat.Y .* (mat.Y - mat.JZ)) * mat.jZ / 2;
Mot=[ Mot' mat.pnul' ];
end

function [ Mot ] = Motifs_pos36(mat)
Mot=(mat.NT .* (mat.MT * (mat.Z.*mat.X))) * mat.jZ;
Mot=[ mat.znul'  Mot' ];
end


function [ Mot ] = Motifs_pos37(mat)
Mot=(mat.MT .* (mat.MT * (mat.Y .* mat.X))) * mat.jZ / 2;
Mot=[ mat.znul'  Mot' ];
end

function [ Mot ] = Motifs_pos38(mat)
Mot=(mat.Z .* mat.Y .* mat.X) * mat.jZ;
Mot=[Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos39(mat)
Mot=(mat.NT .* (mat.MT * (mat.Z .* (mat.Z - mat.JZ)))) * mat.jZ / 2;
Mot=[mat.znul' Mot' ];
end

function [ Mot ] = Motifs_pos40(mat)
Mot=(mat.MT .* (mat.MT * (mat.Y .* (mat.Z - mat.JZ)))) * mat.jZ;
Mot=[mat.znul' Mot' ];
end

function [ Mot ] = Motifs_pos41(mat)
Mot=(mat.Z .* (mat.Z - mat.JZ) .* mat.X) * mat.jZ / 2;
Mot=[ Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos42(mat)
Mot=(mat.Z .* (mat.Z - mat.JZ) .* mat.Y) * mat.jZ / 2;
Mot=[ Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos43(mat)
Mot=(mat.MT .* (mat.MT * ((mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ)))) * mat.jZ / 4 -  mat.MT * ((mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ)) / 4;
Mot=[ mat.znul' Mot' ];
end

function [ Mot ] = Motifs_pos44(mat)
Mot=((mat.Z .* (mat.Z - mat.JZ) .* (mat.Z - 2 * mat.JZ)) * mat.jZ  -  mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 * mat.jZ) )/ 6;
Mot=[  Mot' mat.pnul'];
end

function [ Mot ] = Motifs_pos45(mat)
Mot=mat.MT * ((mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ)) / 6;
Mot=[   mat.znul' Mot'];
end

function [ Mot ] = Motifs_pos46(mat)
Mot=mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ) / 24;
Mot=[ Mot'    mat.pnul'];
end


%%

function [ Mot ] = Motifs_pos47(mat)
Mot=mat.dP.*(mat.dP - mat.jP).*(mat.dP - 2.*mat.jP).*(mat.dP - 3.*mat.jP).*(mat.dP - 4.*mat.jP) / 120;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos48(mat)
Mot=mat.M * ((mat.dP - mat.jP).*(mat.dP - 2.*mat.jP).*(mat.dP - 3.*mat.jP).*(mat.dP - 4.*mat.jP)) / 24;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos49(mat)
Mot=(mat.P.*mat.R.*(mat.R - mat.JP).*(mat.R - 2.*mat.JP)) * mat.jP / 6;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos50(mat)
Mot=(mat.P.*mat.Q.*(mat.Q - mat.JP).*(mat.Q - 2.*mat.JP)) * mat.jP / 6;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos51(mat)
Mot=(mat.N.*(mat.M * (mat.P.*(mat.Q - mat.JP).*(mat.Q - 2.*mat.JP)))) * mat.jP / 2;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos52(mat)
Mot=(mat.M.*(mat.M * (mat.Q.*(mat.Q - mat.JP).*(mat.Q - 2.*mat.JP)))) * mat.jP / 6;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos53(mat)
Mot=(mat.P.*mat.Q.*mat.R.*(mat.R - mat.JP)) * mat.jP / 2;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos54(mat)
Mot=(mat.P.*mat.Q.*(mat.Q - mat.JP).*mat.R) * mat.jP / 2;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos55(mat)
Mot=(mat.M.*(mat.N * (mat.P.*mat.Q.*(mat.Q - mat.JP)))) * mat.jP / 2;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos56(mat)
Mot=(mat.M.*(mat.N * (mat.P.*mat.Q.*(mat.R - mat.JP)))) * mat.jP;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos57(mat)
Mot=(mat.M.*(mat.M * (mat.Q.*(mat.Q - mat.JP).*mat.R))) * mat.jP / 2;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos58(mat)
Mot=(mat.P.*(mat.P - mat.JP).*mat.R.*(mat.R - mat.JP)) * mat.jP / 4;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos59(mat)
Mot=(mat.P.*(mat.P - mat.JP).*mat.Q.*(mat.Q - mat.JP)) * mat.jP / 4;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos60(mat)
Mot=(mat.N.*(mat.M * ((mat.Q - mat.JP).*mat.P.*(mat.P - mat.JP)))) * mat.jP / 2;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos61(mat)
Mot=(mat.M.*(mat.M * ((mat.P - mat.JP).*mat.Q.*(mat.Q - mat.JP)))) * mat.jP / 2;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos62(mat)
Mot=(mat.P.*(mat.P - mat.JP).*mat.Q.*mat.R) * mat.jP / 2;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos63(mat)
Mot=(mat.N.*(mat.M * (mat.P.*(mat.P - mat.JP).*mat.R))) * mat.jP / 2;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos64(mat)
Mot=(mat.M.*(mat.M * ((mat.P - mat.JP).*mat.Q.*mat.R))) * mat.jP / 2;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos65(mat)
Mot=(mat.P.*(mat.P - mat.JP).*(mat.P - 2.*mat.JP).*mat.R) * mat.jP / 6;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos66(mat)
Mot=(mat.P.*(mat.P - mat.JP).*(mat.P - 2.*mat.JP).*mat.Q) * mat.jP / 6;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos67(mat)
Mot=(mat.N.*(mat.M * (mat.P.*(mat.P - mat.JP).*(mat.P - 2.*mat.JP)))) * mat.jP / 6;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos68(mat)
Mot=(mat.M.*(mat.M * ((mat.P - mat.JP).*(mat.P - 2.*mat.JP).*mat.Q))) * mat.jP / 2;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos69(mat)
Mot=(mat.P.*(mat.P - mat.JP).*(mat.P - 2.*mat.JP).*(mat.P - 3.*mat.JP)) * mat.jP / 24 - mat.dP.*(mat.dP - mat.jP).*(mat.dP - 2.*mat.jP).*(mat.dP - 3.*mat.jP) / 24;
Mot=[mat.znul'  Mot'];
end

function [ Mot ] = Motifs_pos70(mat)
Mot=(mat.M.*(mat.M * ((mat.P - mat.JP).*(mat.P - 2.*mat.JP).*(mat.P - 3.*mat.JP)))) * mat.jP / 12 - mat.M * ((mat.dP - mat.jP).*(mat.dP - 2.*mat.jP).*(mat.dP - 3.*mat.jP)) / 12;
Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos71(mat)
Mot=sum(sum((mat.NTB .* (mat.NTB - mat.JP3) .* mat.MTA .* mat.KP3),2),3)/2;
Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos72(mat)
  Mot=sum(sum((mat.MTD .* (mat.MTD - mat.JP3) .* mat.MTA .* mat.KP3),2),3) / 4;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos73(mat)
  Mot=sum(sum((mat.NB .* (mat.NB - mat.JZ3) .* mat.MA .* mat.KZ3),2),3) / 2;
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos74(mat)
  Mot=sum(sum((mat.MD .* (mat.MD - mat.JZ3) .* mat.MA .* mat.KZ3),2),3) / 4;
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos75(mat)
  Mot=sum(sum(mat.MTD .* mat.MTA .* mat.NTB .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos76(mat)
  Mot=sum(sum(mat.MTA .* mat.NTB .* mat.NTC .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos77(mat)
  Mot=sum(sum(mat.MD .* mat.MB .* mat.MC .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos78(mat)
  Mot=sum(sum(mat.MB .* mat.NB .* mat.Na .* mat.KZ3,2),3);
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos79(mat)
  Mot=sum(sum(mat.MTB .* mat.MTC .* mat.MTD .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos80(mat)
  Mot=sum(sum(mat.MTB .* mat.NTB .* mat.NTA .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos81(mat)
  Mot=sum(sum(mat.MA .* mat.MD .* mat.NB .* mat.KZ3,2),3);
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos82(mat)
  Mot=sum(sum(mat.MA .* mat.NB .* mat.NC .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos83(mat)
  Mot=sum(sum(mat.MTB .* mat.NTA .* mat.NTC .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos84(mat)
  Mot=sum(sum(mat.NTA .* mat.MTB .* mat.MTD .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos85(mat)
  Mot=sum(sum(mat.MTB .* mat.MTC .* mat.NTC .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos86(mat)
  Mot=sum(sum(mat.MB .* mat.Na .* mat.NC .* mat.KZ3,2),3);
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos87(mat)
  Mot=sum(sum(mat.MD .* mat.MB .* mat.Na .* mat.KZ3,2),3);
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos88(mat)
  Mot=sum(sum(mat.NB .* mat.MB .* mat.MC .* mat.KZ3,2),3);
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos89(mat)
  Mot=sum(sum(mat.MTA .* mat.NTA .* mat.NTB .* mat.KP3,2),3);
 Mot=[  mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos90(mat)
  Mot=sum(sum(mat.MTA .* mat.MTB .* mat.NTB .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos91(mat)
  Mot=sum(sum(mat.MTA .* mat.MTB .* mat.MTD .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos92(mat)
  Mot=sum(sum(mat.MA .* mat.Na .* mat.NB .* mat.KZ3,2),3);
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos93(mat)
  Mot=sum(sum(mat.MB .* mat.NB .* mat.MA .* mat.KZ3,2),3);
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos94(mat)
  Mot=sum(sum(mat.MA .* mat.MB .* mat.MD .* mat.KZ3,2),3);
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos95(mat)
  Mot=sum(sum(mat.MTB .* mat.NTA .* (mat.NTA - mat.JP3) .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos96(mat)
  Mot=sum(sum(mat.MTB .* (mat.MTB - mat.JP3) .* mat.NTA .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos97(mat)
  Mot=sum(sum(mat.MTB .* mat.MTC .* (mat.MTC - mat.JP3) .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos98(mat)
  Mot=sum(sum(mat.MD .* mat.MA .* mat.Na .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos99(mat)
  Mot=sum(sum(mat.MA .* mat.NB .* mat.MC .* mat.KZ3,2),3);
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos100(mat)
  Mot=sum(sum(mat.MTA .* mat.NTA .* (mat.NTA - mat.JP3) .* mat.KP3,2),3) / 4;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos101(mat)
  Mot=sum(sum(mat.MTA .* mat.MTB .* (mat.MTB - mat.JP3) .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos102(mat)
  Mot=sum(sum(mat.MA .* (mat.MA - mat.JZ3) .* mat.NB .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos103(mat)
  Mot=sum(sum(mat.MD .* mat.MA .* (mat.MA - mat.JZ3) .* mat.KZ3,2),3) / 4;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos104(mat)
  Mot=sum(sum(mat.MTD .* mat.MTA .* mat.NTA .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos105(mat)
  Mot=sum(sum(mat.MTB .* mat.MTA .* mat.NTC .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos106(mat)
  Mot=sum(sum(mat.MB .* mat.Na .* (mat.Na - mat.JZ3) .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos107(mat)
  Mot=sum(sum(mat.MB .* (mat.MB - mat.JZ3) .* mat.Na .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos108(mat)
  Mot=sum(sum(mat.MB .* (mat.MB - mat.JZ3) .* mat.MC .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos109(mat)
  Mot=sum(sum(mat.MTA .* (mat.MTA - mat.JP3) .* mat.NTB .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos110(mat)
  Mot=sum(sum(mat.MTA .* (mat.MTA - mat.JP3) .* mat.MTD .* mat.KP3,2),3) / 4;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos111(mat)
  Mot=sum(sum(mat.MA .* mat.Na .* (mat.Na - mat.JZ3) .* mat.KZ3,2),3) / 4;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos112(mat)
  Mot=sum(sum(mat.MA .* mat.MB .* (mat.MB - mat.JZ3) .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos113(mat)
  Mot=sum(sum(mat.MTB .* mat.MTC .* mat.NTA .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos114(mat)
  Mot=sum(sum(mat.MB .* mat.MC .* mat.Na .* mat.KZ3,2),3) / 2;
 Mot=[Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos115(mat)
  Mot=sum(sum(mat.MTA .* mat.MTB .* mat.NTA .* mat.KP3,2),3);
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos116(mat)
  Mot=sum(sum(mat.MTA .* mat.MTB .* mat.MTC .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos117(mat)
  Mot=sum(sum(mat.MA .* mat.Na .* mat.MB .* mat.KZ3,2),3);
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos118(mat)
  Mot=sum(sum(mat.MA .* mat.MB .* mat.MC .* mat.KZ3,2),3) / 2;
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos119(mat)
  Mot=sum(sum(mat.MTA .* (mat.MTA - mat.JP3) .* mat.NTA .* mat.KP3,2),3) / 4;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos120(mat)
  Mot=sum(sum(mat.MTA .* (mat.MTA - mat.JP3) .* mat.MTB .* mat.KP3,2),3) / 2;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos121(mat)
  Mot=sum(sum(mat.MA .* (mat.MA - mat.JZ3) .* mat.Na .* mat.KZ3,2),3) / 4;
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos122(mat)
  Mot=sum(sum(mat.MB .* mat.MA .* (mat.MA - mat.JZ3) .* mat.KZ3,2),3) / 2;
 Mot=[ Mot'    mat.pnul'];
end

function [ Mot ] = Motifs_pos123(mat)
  Mot=sum(sum(mat.MTA .* (mat.MTA - mat.JP3) .* (mat.MTA - 2 .* mat.JP3) .* mat.KP3,2),3) / 12;
 Mot=[ mat.znul'    Mot'];
end

function [ Mot ] = Motifs_pos124(mat)
  Mot=sum(sum(mat.MA .* (mat.MA - mat.JZ3) .* (mat.MA - 2 .* mat.JZ3) .* mat.KZ3,2),3) / 12;
 Mot=[ Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos125(mat)
  Mot=(mat.NT .* (mat.MT * (mat.Z .* (mat.Y - mat.JZ) .* (mat.Y - 2 .* mat.JZ)))) * mat.jZ / 2;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos126(mat)
  Mot=(mat.MT .* (mat.MT * (mat.Y .* (mat.Y - mat.JZ) .* (mat.Y - 2 .* mat.JZ)))) * mat.jZ / 6;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos127(mat)
  Mot=(mat.Z .* mat.X .* (mat.X - mat.JZ) .* (mat.X - 2 .* mat.JZ)) * mat.jZ / 6;
 Mot=[ Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos128(mat)
  Mot=(mat.Z .* mat.Y .* (mat.Y - mat.JZ) .* (mat.Y - 2 .* mat.JZ)) * mat.jZ / 6;
 Mot=[ Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos129(mat)
  Mot=(mat.MT .* (mat.NT * (mat.Z .* mat.Y .* (mat.Y - mat.JZ)))) * mat.jZ / 2;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos130(mat)
  Mot=(mat.MT .* (mat.NT * (mat.Z .* mat.Y .* (mat.X - mat.JZ)))) * mat.jZ;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos131(mat)
  Mot=(mat.MT .* (mat.MT * (mat.Y .* (mat.Y - mat.JZ) .* mat.X))) * mat.jZ / 2;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos132(mat)
  Mot=(mat.Z .* mat.Y .* mat.X .* (mat.X - mat.JZ)) * mat.jZ / 2;
 Mot=[ Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos133(mat)
  Mot=(mat.Z .* mat.Y .* (mat.Y - mat.JZ) .* mat.X) * mat.jZ / 2;
 Mot=[  Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos134(mat)
  Mot=(mat.NT .* (mat.MT * ((mat.Y - mat.JZ) .* mat.Z .* (mat.Z - mat.JZ)))) * mat.jZ / 2;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos135(mat)
  Mot=(mat.MT .* (mat.MT * ((mat.Z - mat.JZ) .* mat.Y .* (mat.Y - mat.JZ)))) * mat.jZ / 2;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos136(mat)
  Mot=(mat.Z .* (mat.Z - mat.JZ) .* mat.X .* (mat.X - mat.JZ)) * mat.jZ / 4;
 Mot=[  Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos137(mat)
  Mot=(mat.Z .* (mat.Z - mat.JZ) .* mat.Y .* (mat.Y - mat.JZ)) * mat.jZ / 4;
 Mot=[  Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos138(mat)
  Mot=(mat.NT .* (mat.MT * (mat.Z .* (mat.Z - mat.JZ) .* mat.X))) * mat.jZ / 2;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos139(mat)
  Mot=(mat.MT .* (mat.MT * ((mat.Z - mat.JZ) .* mat.Y .* mat.X))) * mat.jZ / 2;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos140(mat)
  Mot=(mat.Z .* (mat.Z - mat.JZ) .* mat.Y .* mat.X) * mat.jZ / 2;
 Mot=[  Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos141(mat)
  Mot=(mat.NT .* (mat.MT * (mat.Z .* (mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ)))) * mat.jZ / 6;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos142(mat)
  Mot=(mat.MT .* (mat.MT * ((mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ) .* mat.Y))) * mat.jZ / 2;
 Mot=[ mat.znul'    Mot'];
end


  function [ Mot ] = Motifs_pos143(mat)
  Mot=(mat.Z .* (mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ) .* mat.X) * mat.jZ / 6;
 Mot=[ Mot'    mat.pnul'];
end


  function [ Mot ] = Motifs_pos144(mat)
  Mot=(mat.Z .* (mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ) .* mat.Y) * mat.jZ / 6;
 Mot=[  Mot'    mat.pnul'];
  end

  function [ Mot ] = Motifs_pos145(mat)
  Mot=(mat.MT .* (mat.MT * ((mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ) .* (mat.Z - 3 .* mat.JZ)))) * mat.jZ / 12 - mat.MT * ((mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ)) / 12;
 Mot=[ mat.znul'    Mot'];
  end

 function [ Mot ] = Motifs_pos146(mat)
  Mot=(mat.Z .* (mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ) .* (mat.Z - 3 .* mat.JZ)) * mat.jZ / 24 - mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ) / 24;
 Mot=[ Mot'    mat.pnul'];
end

  function [ Mot ] = Motifs_pos147(mat)
  Mot=mat.MT * ((mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ) .* (mat.dZ - 4 .* mat.jZ)) / 24;
 Mot=[ mat.znul'    Mot'];
  end

  function [ Mot ] = Motifs_pos148(mat)
  Mot=mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ) .* (mat.dZ - 4 .* mat.jZ) / 120;
 Mot=[ Mot'    mat.pnul'];
end