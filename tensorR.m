function [ Rx ] = tensorR( A, B)

dimA=size(A);
dimB=size(B);

alongA=2;
alongB=1;


seqA=1:dimA;
AllA=false;
seqAx=seqA;
seqAx(alongA)=[];
permA=[seqAx,alongA];
dimAx=dimA;
dimAx(alongA)=[];
Ax=reshape(A,[prod(dimAx),prod(dimA(alongA))]);


seqB=1:dimB;
AllB=false;
seqBx=seqB;
seqBx(alongB)=[];
permB=[alongB,seqBx];
dimBx=dimB;
dimBx(alongB)=[];
Bx=reshape(B,[prod(dimB(alongB)),prod(dimBx)]);


R=Ax*Bx;

Rx=reshape(R,[dimAx,dimBx]);

end