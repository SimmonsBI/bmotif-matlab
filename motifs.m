function [ Motifout ] = Motifs(M,V_motifs)
% Counts the frequency with which bipartite motifs up to six nodes occur 
% in a network
%   M is a binary two-dimensional matrix
%   V_motifs: is a vector with the id of motifs you want count (between 1
%   and 44)
% Output: a table with one row for each motif and two columns. The first column 
% (IDmotif) gives the ID of the motif, the second column gives the frequency with which 
% that motif occurs in the network

if nargin==1
    if Check_mot(M)==0
        return
    end
elseif nargin==2
    if Check_mot(M,V_motifs)==0
        return
    else
        Mot_var_id={'ID','Motifs'};
    end
elseif nargin>2
    print 'insert the correct number of arguments'
    return
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


mat.MTAK = mat.MTA.*mat.KP3;
mat.MTBK = mat.MTB.*mat.KP3;
mat.MTCK = mat.MTC.*mat.KP3;
mat.MTDK = mat.MTD.*mat.KP3;

mat.MAK = mat.MA.*mat.KZ3;
mat.MBK = mat.MB.*mat.KZ3;
mat.MCK = mat.MC.*mat.KZ3;
mat.MDK = mat.MD.*mat.KZ3;

mat.NTAK = mat.NTA.*mat.KP3;
mat.NTBK = mat.NTB.*mat.KP3;
mat.NTCK = mat.NTC.*mat.KP3;

mat.NaK = mat.Na.*mat.KZ3;
mat.NBK = mat.NB.*mat.KZ3;
mat.NCK = mat.NC.*mat.KZ3;




%Mot_name={'ID','#Motids'};
count=1;
for n_motifs=V_motifs

    Motif=0;
    %%Motifs 1 count
    
    fh=str2func(strcat('Motifs',num2str(n_motifs)));
    Motif=fh(mat);
    
    Motifout(count,:)=[n_motifs,Motif];
    count=count+1;
end

Motifout=array2table(Motifout,'VariableNames',Mot_var_id);

%Motifout=array2table(Motifv)
end


%%
%%MOTIFS Functions
function [ Mot ] = Motifs1(mat)
Mot=sum(sum(mat.M));
end


function [ Mot ] = Motifs2(mat)
Mot= sum(mat.dZ .* (mat.dZ - mat.jZ)) / 2;

end

function [ Mot ] = Motifs3(mat)
Mot= sum(mat.dP .* (mat.dP - mat.jP)) / 2;
end


function [ Mot ] = Motifs4(mat)
Mot= sum(mat.dP .* (mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP)) / 6;
end

function [ Mot ] = Motifs5(mat)
Mot=sum(sum(mat.Z .* mat.Y));
end

function [ Mot ] = Motifs6(mat)
Mot=sum(sum(mat.Z .* (mat.Z - mat.JZ))) / 4 - sum(sum(mat.dZ .* (mat.dZ - mat.jZ))) / 4;
end

function [ Mot ] = Motifs7(mat)
Mot= sum(sum(mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ))) / 6;
end

function [ Mot ] = Motifs8(mat)
Mot=sum(sum(mat.dP .* (mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP) .* (mat.dP - 3 .* mat.jP))) / 24;
end

function [ Mot ] = Motifs9(mat)
Mot= sum(sum(mat.P .* mat.Q .* (mat.Q - mat.JP))) / 2;
end

function [ Mot ] = Motifs10(mat)
Mot=sum(sum(mat.P .* mat.Q .* mat.R)) / 2;
end

function [ Mot ] = Motifs11(mat)
Mot=sum(sum(mat.P .* (mat.P - mat.JP) .* mat.Q)) / 2;
end

function [ Mot ] = Motifs12(mat)
Mot=sum(sum(mat.P .* (mat.P - mat.JP) .* (mat.P - 2 .* mat.JP))) / 12 - sum(sum(mat.dP .* (mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP))) / 12;
end

function [ Mot ] = Motifs13(mat)
Mot=sum(sum(mat.Z .* mat.Y .* (mat.Y - mat.JZ))) / 2;
end

function [ Mot ] = Motifs14(mat)
Mot=sum(sum(mat.Z .* mat.Y .* mat.X)) / 2;
end

function [ Mot ] = Motifs15(mat)
Mot=sum(sum(mat.Z .* (mat.Z - mat.JZ) .* mat.Y)) / 2;
end

function [ Mot ] = Motifs16(mat)
Mot= sum(sum(mat.Z .* (mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ))) / 12 - sum(sum(mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ))) / 12;
end

function [ Mot ] = Motifs17(mat)
Mot= sum(sum(mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ))) / 24;
end

function [ Mot ] = Motifs18(mat)
Mot= sum(sum(mat.dP .* (mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP) .* (mat.dP - 3 .* mat.jP) .* (mat.dP - 4 .* mat.jP)) ) / 120;
end

function [ Mot ] = Motifs19(mat)
Mot= sum(sum(mat.P .* mat.Q .* (mat.Q - mat.JP) .* (mat.Q - 2 .* mat.JP)) ) / 6;
end

function [ Mot ] = Motifs20(mat)
Mot= sum(sum(mat.Q .* mat.P .* mat.R .* (mat.Q - mat.JP)) ) / 2;
end

function [ Mot ] = Motifs21(mat)
Mot= sum(sum(mat.Q .* (mat.Q - mat.JP) .* mat.P .* (mat.P - mat.JP)) ) / 4;
end

function [ Mot ] = Motifs22(mat)
Mot= sum(sum(mat.Q .* mat.P .* (mat.P - mat.JP) .* mat.R) ) / 4;
end

function [ Mot ] = Motifs23(mat)
Mot= sum(sum(mat.P .* (mat.P - mat.JP) .* (mat.P - 2 .* mat.JP) .* mat.Q) ) / 6;
end

function [ Mot ] = Motifs24(mat)
Mot= sum(sum(mat.P .* (mat.P - mat.JP) .* (mat.P - 2 .* mat.JP) .* (mat.P - 3 .* mat.JP))) / 48 - sum(sum(mat.dP .* (mat.dP - mat.jP) .* (mat.dP - 2 .* mat.jP) .* (mat.dP - 3 .* mat.jP))) / 48;
end

function [ Mot ] = Motifs25(mat)
Mot= sum(sum(sum(mat.MAK .* mat.NBK .* (mat.NBK - mat.JZ3))) ) / 4;
end

function [ Mot ] = Motifs26(mat)
Mot=  sum(sum(sum(mat.MDK .* mat.MBK .* mat.MCK))) / 2;
end

function [ Mot ] = Motifs27(mat)
Mot=  sum(sum(sum(mat.NBK .* mat.NCK .* mat.MAK))) / 2;
end

function [ Mot ] = Motifs28(mat)
Mot=  sum(sum(sum(mat.MBK .* mat.MCK .* mat.NCK)));
end

function [ Mot ] = Motifs29(mat)
Mot=  sum(sum(sum(mat.MAK .* mat.MBK .* mat.MDK)));
end


function [ Mot ] = Motifs30(mat)
Mot=  sum(sum(sum(mat.MAK .* mat.MBK .* mat.NCK)))/2;
end

function [ Mot ] = Motifs31(mat)
Mot=  sum(sum(sum(mat.MAK .* mat.NBK .*(mat.MAK - mat.JZ3))))/4;
end

function [ Mot ] = Motifs32(mat)
Mot=  sum(sum(sum(mat.MBK .*( mat.MBK - mat.JZ3).*mat.MCK)))/2;
end

function [ Mot ] = Motifs33(mat)
Mot=  sum(sum(sum(mat.MBK .*( mat.MBK - mat.JZ3).*mat.MAK)))/4;
end

function [ Mot ] = Motifs34(mat)
Mot=  sum(sum(sum(mat.NaK .* mat.MBK .* mat.MCK)))/6;
end

function [ Mot ] = Motifs35(mat)
Mot=  sum(sum(sum(mat.MAK .* mat.MBK .* mat.MCK)))/2;
end

function [ Mot ] = Motifs36(mat)
Mot=  sum(sum(sum(mat.MAK .*( mat.MAK -mat.JZ3) .* mat.MBK)))/4;
end

function [ Mot ] = Motifs37(mat)
Mot=  sum(sum(sum(mat.MAK .*( mat.MAK -mat.JZ3) .* (mat.MAK -2.*mat.JZ3))))/36;
end

function [ Mot ] = Motifs38(mat)
Mot=  sum(sum(sum(mat.Z .*mat.Y.*( mat.Y -mat.JZ) .* (mat.Y-2.*mat.JZ))))/6;
end

function [ Mot ] = Motifs39(mat)
Mot=  sum(sum(sum(mat.Z .*( mat.Y -mat.JZ) .* mat.X .*mat.Y)))/2;
end

function [ Mot ] = Motifs40(mat)
Mot=  sum(sum(sum(mat.Z .*( mat.Z -mat.JZ) .* mat.Y .*(mat.Y-mat.JZ))))/4;
end

function [ Mot ] = Motifs41(mat)
Mot=  sum(sum(sum(mat.Z .*( mat.Z -mat.JZ) .* mat.X .* mat.Y)))/4;
end

function [ Mot ] = Motifs42(mat)
Mot=  sum(sum(sum(mat.Z .*( mat.Z -mat.JZ).* (mat.Z - 2 .* mat.JZ) .* mat.Y)))/6;
end

function [ Mot ] = Motifs43(mat)
Mot= sum(sum(sum(mat.Z .* (mat.Z - mat.JZ) .* (mat.Z - 2 .* mat.JZ) .* (mat.Z - 3 .* mat.JZ)))) / 48 -   sum(sum(sum(mat.dZ .* (mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ)))) / 48;
end

function [ Mot ] = Motifs44(mat)
Mot=  sum(sum(sum(mat.dZ .*( mat.dZ - mat.jZ) .* (mat.dZ - 2 .* mat.jZ) .* (mat.dZ - 3 .* mat.jZ) .* (mat.dZ - 4 .* mat.jZ))))/120;
end

