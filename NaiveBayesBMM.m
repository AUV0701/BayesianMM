
%%% SETS THE VALUE OF THETA* %%%
p = randi([1,30],2,1);
p = p/sum(p);
%%% SETS RANDOMLY THE VALUES OF A, B, C %%%
A = randi([10,30],3,1);
A = A/sum(A);
A = [1,0,0];
%%%%%%%%%%%%%%%%%%%
%%% NAIVE BAYES %%%
%%%%%%%%%%%%%%%%%%%

%%% SET THE VALUES OF phi_1 and phi_2 %%%
P(:,1) = [0.01 : 0.01 : 1];
phi2(:,1) = [0.01 : 0.01 : 1];
for i = 1: size(P,1)-1
    phi2 = vertcat(phi2,circshift(P,i));
end
phi1(:,1) = repmat(P,size(P,1),1);
%%% COMPUTES (1 - phi_1) AND (1 - phi_2) %%%
PHI1 = 1 - phi1;
PHI2 = 1 - phi2;

%%% CALCULATE THE INDIVIDUAL TERMS FOR THE EXPECTED VALUE FOR EACH SET OF phi_1 and phi_2%%%
term1 = (p(1)*phi1 + p(2)*phi2).*phi1.*phi2./(A(1)*phi1.^2 + A(2)*phi2.^2 + A(3)*phi1.*phi2);
term2 = (p(1)*PHI1 + p(2)*PHI2).*PHI1.*PHI2./(A(1)*PHI1.^2 + A(2)*PHI2.^2 + A(3)*PHI1.*PHI2);
%%% CACLULATE THE EXPECTED VALUE AND PLOT %%%
delta = 1 - (term1 + term2);
DELTA = vec2mat(delta,size(P,1));
x = phi1(1:size(P,1),1);
y = phi2(1:size(P,1),1);
figure(2);surf(x,y,DELTA')
xlabel('\phi_{1}');
ylabel('\phi_{2}')

% %%
% phi = [0.001 : 0.001 : 0.02];
% k = [10:1:700];
% for i = 1 : size(phi,2)
%     for j = 1 : size(k,2)
%         term1(i,j) = k(j)*phi(i)*(p(1)*k(j) + p(2))/( A(1)*k(j)^2 + A(2) + A(3)*k(j) );
%         term2(i,j) = (1 - k(j)*phi(i))*(1-phi(i))*(1 - phi(i)*(p(1)*k(j) + p(2)))/( 1 - phi(i)^2*(A(1)*k(j)^2 + A(2) - k(j)*A(3)) - phi(i)*A(3)*(1 -k(j)) );
%     end
% end
% 
% Delta = 1 - (term1 + term2);
% surf(phi,k,Delta')
% figure(2);surf(phi,k,Delta')

