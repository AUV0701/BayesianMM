function [dataStream, p, mu, sigma, probability] = generateData(numTrials,nComp)

% WEIGHTS OF THE GAUSSIAN COMPONENTS
p = randi([3,10],nComp,1);
p = p/sum(p);
% PARAMETERS TO GENERATE DATA
mu = [2, 3];
sigma = [1, 1];
for i = 1 : nComp
    cases = ceil(p(i)*numTrials);
    data{i} = normrnd(mu(i), sigma(i), cases,1);
end
dataStream = vertcat(data{1}, data{2});
dataStream = dataStream(randperm(size(dataStream,1)));
probability = p(1)*mvnpdf(dataStream, mu(1),sigma(1)) + p(2)*mvnpdf(dataStream, mu(2),sigma(2));
% %%
% subplot(1,3,1)
% axis square
% s = scatter(accuracyBL1',accuracyFinal,sz,linspace(20,100,length(accuracyFinal)),'filled')
% hold on; plot([1:7:100], [1:7:100], '--k')
% 
% subplot(1,3,2)
% axis square
% s = scatter(accuracyEM',accuracyFinal,sz,linspace(20,100,length(accuracyFinal)),'filled')
% hold on; plot([1:7:100], [1:7:100], '--k')
% 
% subplot(1,3,3)
% axis square
% s = scatter(SleepStageResults,accuracyFinal,sz,linspace(20,100,length(accuracyFinal)),'filled')
% hold on; plot([1:7:100], [1:7:100], '--k')
