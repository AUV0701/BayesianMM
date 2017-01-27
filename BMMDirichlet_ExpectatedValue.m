
p = randi([3,100],3,1);
p = p/sum(p);

c = randi([1,20],2,1);
c = c/sum(c);

mu = [2,2.2];
sigma = [1,1];

x = rand(1000000,1);
x = bsxfun(@times, 30, x) -15;

pdf(:,1) = mvnpdf(dataStream,mu(1),sigma(1));
pdf(:,2) = mvnpdf(dataStream,mu(2),sigma(2));

f = 1 - (pdf(:,1).*pdf(:,2)).*(c(1)*pdf(:,1) + c(2)*pdf(:,2))./(p(1)*pdf(:,1).^2 + p(2)*pdf(:,2).^2 + p(3)*pdf(:,1).*pdf(:,2));

plot(x,f,'*');

%%
alpha = abs(rand(1,nComp)) + 1;
A = zeros(size(dataStream,1),1);
B = zeros(size(dataStream,1),1);
ALPHA = zeros(size(dataStream,1),2);
for i = 1 : size(dataStream,1)
    ALPHA(i,:) = alpha;
    alpha(1,1) = alpha(1,1) + (alpha(1,1)*(alpha(1,1)+1)*pdf(i,1)*(pdf(i,1) - pdf(i,2)))/ (alpha(1,1)*(alpha(1,1)+1)*(pdf(i,1))^2 + alpha(1,2)*(alpha(1,2)+1)*(pdf(i,2))^2 +2*((alpha(1,2)+1)*(alpha(1,1)+1)*pdf(i,1)*pdf(i,2))) ;
    alpha(1,2) = alpha(1,2) + (alpha(1,2)*(alpha(1,2)+1)*pdf(i,2)*(pdf(i,2) - pdf(i,1)))/ (alpha(1,1)*(alpha(1,1)+1)*(pdf(i,1))^2 + alpha(1,2)*(alpha(1,2)+1)*(pdf(i,2))^2 +2*((alpha(1,2)+1)*(alpha(1,1)+1)*pdf(i,1)*pdf(i,2))) ;
    
    A(i) =  alpha(1,1)*( alpha(1,1) + 1)/((sum(alpha)+1)*(sum(alpha)+2));
    B(i) =  alpha(1,2)*( alpha(1,2) + 1)/((sum(alpha)+1)*(sum(alpha)+2));
end