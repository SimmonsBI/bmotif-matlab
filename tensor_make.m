function [ output ] = tensor_make( Ma,Mb )
x=size(Ma,1);
y=size(Ma,2);
z=size(Mb,2);

output=ones(x,y,z);
for j=1:y
    for k=1:z
        output(:,j,k)=Ma(:,j).*Mb(:,k);
    end
end
end
