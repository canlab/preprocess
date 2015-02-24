function [M] = latin_square(n)


if n>30
    warning('This might take a minute...');
elseif n>50
    warning('This will probably take a while...');
elseif n>80
    warning('This is going to take a really long time...');
end

M=zeros(n);

M(1,:)=randperm(n);

i=2;
while i<=n
    for j=1:n
        a=M(1:i-1,j);
        b=M(i,1:j-1);
        for k=randperm(n)
            if ~sum(a==k)&&~sum(b==k)
                M(i,j)=k;
                break
            end
        end
        if ~M(i,j)
            M(i,:)=0;
            i=i-1;
            break
        end
    end
    i=i+1;
end