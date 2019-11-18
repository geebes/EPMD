clf

n=10;

kern=diag(0.9.*ones(1,n),0)+diag(0.05.*ones(1,n-1),1)+diag(0.05.*ones(1,n-1),-1);


colormap(flipud(hot))

subplot(221)
x=zeros(n)./n;
x(5,5)=n;
for j=1:5
    x=kern*x*kern;
end
pcolor(padarray(x,[1 1],'post'))
set(gca,'XTick',1+(0.5:n),'XTickLabel',1:n,...
        'YTick',1+(0.5:n),'YTickLabel',1:n)
xlabel('Phenotype')
ylabel('Location')

subplot(222)
for i=1:n
    x(:,i)=normpdf(1:n,i,1.5);
end
pcolor(padarray(x,[1 1],'post'))
set(gca,'XTick',1+(0.5:n),'XTickLabel',1:n,...
        'YTick',1+(0.5:n),'YTickLabel',1:n)
xlabel('Phenotype')
ylabel('Location')

subplot(223)
xnew=sparse(zeros(n,n^2));
for i=1:n % for each phenotype
    xnew(:,i+(i-1).*n)=x(:,i);
    YTck{i}=[num2str(i) '.' num2str(i)];
end
x=xnew;
pcolor(padarray(x,[1 1],'post'))
set(gca,'XTick',(1:n)+((1:n)-1).*n,'XTickLabel',YTck,...
        'YTick',1+(0.5:n),'YTickLabel',1:n)
xlabel('Phenotype')
ylabel('Location')
    
subplot(224)
for i=1:n
    for j=1:10
        x(:,(1:n)+(i-1).*n)=x(:,(1:n)+(i-1).*n)*kern;
    end
end
pcolor(padarray(x,[1 1],'post'))
set(gca,'XTick',(1:n)+((1:n)-1).*n,'XTickLabel',YTck,...
        'YTick',1+(0.5:n),'YTickLabel',1:n)
xlabel('Phenotype')
ylabel('Location')









