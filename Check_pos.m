function [ flag_e ] = Check_pos(M,V_motif)
flag_e=1;

if iscell(M)==1
    flag_e=0;
    sprintf 'M is not a matrix'
    return
end

if isnumeric(M)==0
    flag_e=0;
    sprintf 'M is not a matrix'
    return
end

if ismatrix(M)==0
    flag_e=0;
    sprintf 'M is not a matrix'
    return
end

if size(M,1)==1 & size(M,2)==1
    flag_e=0;
    sprintf 'M has only one value'
    return
end


if length(find(M==0))+length(find(M==1))~=size(M,1)*size(M,2)
    flag_e=0;
    sprintf 'M has values different from 1 and 0'
    return
end

z=size(M,1);
p=size(M,2);

if nargin<2
    sprintf 'missing which motif do you want count'
    flag_e=0;
    return
    
end

if nargin==2
    if isvector(V_motif)~=1
        flag_e=0;
        sprintf 'motifs requested is not vector'
        return
    else
        for n_motif=V_motif
            if isnumeric(n_motif)==0
                flag_e=0;
                sprintf 'number motifs is not a number'
                return
            end
            if n_motif<1 || n_motif>946
                flag_e=0;
                sprintf 'number of position motifs not valid enter a number between [1:148]'
                return
            end
            if mod(n_motif,1)~=0
                flag_e=0;
                sprintf 'number of motifs not valid enter a integer number between [1:148]'
                return
            end
        end 
    end
end
    
end

