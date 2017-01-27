%% BMM Learning
alpha = abs(rand(1,nComp)) + 1;
%mu = [2, 3];
%sigma = [1, 1];
k = 1;
for i = 1 : size(dataStream)
    hyperParameters(i,:) = alpha;
    expected_hyperParameters(i,:) = mean(hyperParameters,1);
    sumHyperparameters(i) = sum(alpha);
    sum_expected_hyperParameters(i,:) = mean(sumHyperparameters);
    x = dataStream(i);
    for j = 1 : nComp
        p(j) = mvnpdf(x, mu(j), sigma(j));
    end
    normalisationConstant = sum(alpha.*reshape(p,size(alpha)));
    w = (alpha.*reshape(p,size(alpha)))./normalisationConstant; 
   %%% Evaluate Moments %%%
   constant = sum(alpha)+1;
   for j = 1 : nComp
       moment1(i,j) = (w(j)*(alpha(j)+1) + (1-w(j))*alpha(j))/constant;
       moment2(i,j) = ( w(j)*(alpha(j) + 1)*(alpha(j) + 2) + (1-w(j))*alpha(j)*(alpha(j)+1) )/(constant*(constant+1));
   end
   alpha = (moment1(i,:).*(moment1(i,:) - moment2(i,:)))./(moment2(i,:) - moment1(i,:).^2);
   if (sum(alpha) < sumHyperparameters(i))
       dataContradiction(k,1) = i;
       dataContradiction(k,2) = x;
       parametersContradiction(k,:) = hyperParameters(i,:);
       parametersContradictionFinal(k,:) = alpha;
       k = k+1;
   end
end
%% BMM Learning Expected
alpha = abs(rand(1,nComp)) + 1;
%mu = [2, 3];
%sigma = [1, 1];
k = 1;
for i = 1 : size(dataStream)
    hyperParameters(i,:) = alpha;
    sumHyperparameters(i) = sum(alpha);
    x = dataStream(i);
    for j = 1 : nComp
        p(:,j) = mvnpdf(dataStream, mu(j), sigma(j));
    end
    %normalisationConstant = sum(alpha.*reshape(p,size(alpha)));
    normalisationConstant = sum(bsxfun(@times,alpha,p),2);
    %w = (alpha.*reshape(p,size(alpha)))./normalisationConstant; 
    w = bsxfun(@ldivide, normalisationConstant,bsxfun(@times,alpha,p));
   %%% Evaluate Moments %%%
   constant = sum(alpha)+1;
   for j = 1 : nComp
       moment1(:,j) = (bsxfun(@times,(alpha(j)+1), w(:,j)) + bsxfun(@times,alpha(j), (1-w(:,j))))/constant;
       moment2(:,j) = (bsxfun(@times,((alpha(j) + 1)*(alpha(j) + 2)), w(:,j)) + bsxfun(@times,(alpha(j)*(alpha(j)+1)), (1-w(:,j))))/(constant*(constant+1));
   end
   alphaExpected = (moment1(:,:).*(moment1(:,:) - moment2(:,:)))./(moment2(:,:) - moment1(:,:).^2);
   alpha = alphaExpected(i,:);
   expectedHyperparameters(i,:) = mean(alphaExpected,1);
   if (sum(alpha) < sumHyperparameters(i))
       dataContradiction(k,1) = i;
       dataContradiction(k,2) = x;
       parametersContradiction(k,:) = hyperParameters(i,:);
       parametersContradictionFinal(k,:) = alpha;
       k = k+1;
   end
end
%% plot of hyperparameters for each iteration
plot(hyperParameters);
hold on; plot(sumHyperparameters,'k')
axis square
xlim([0,size(dataStream,1)])
xlabel('# Observations')
ylabel('Value')
title('hyper-parameter values')
legend('\alpha_1','\alpha_2','\alpha_1 + \alpha_2')

%% scatter plot of data points
plot(dataStream,'r*');
hold on; plot(dataContradiction(:,1),dataContradiction(:,2),'ks')
axis square
xlim([0,size(dataStream,1)])
xlabel('# Observations')
ylabel('Value')
title('scatter plot of observations')
legend('data','contradictory points')
%% density plot of data points
plot(dataStream, probability,'r*')
hold on; plot(dataContradiction(:,2),probability(dataContradiction(:,1)),'ks')
axis square
xlabel('Observation')
ylabel('density')
title('density plot (\mu_1 = 2 , \mu_2 = 4.5, \sigma_1^2 = \sigma_2^2 = 0.5, w_1 = 0.667, w_2 = 0.333)')
legend('data','contradictory points')