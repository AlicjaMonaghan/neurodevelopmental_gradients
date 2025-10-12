function tab = dominance(X,y)
%% function for conducitng domiance analysis based on multiple regression
%%INPUT:
%X: matrix with columns holding the explanatory variables
%y: vector holding the response variable 
%%OUTPUT
%tb: table sumarizing the results 
%%
% remove missing data
N = size(X,2);
data = rmmissing([X y],1); 
X = data(:,1:N);
y = data(:,N+1);
%calculate R2 for the complete linear model (y = X_1+...+X_N)
mdl = fitlm(X,y);
R2_full = mdl.Rsquared.Ordinary; 
%calculte R2 for individual linear models (y = X_i, i = 1:N)
R2_indi = zeros(N,1);
for i = 1:N
    mdl = fitlm(X(:,i),y);
    R2_indi(i) = mdl.Rsquared.Ordinary;
end
%calculate R2 for marginal linear models (y = X_k+...,+X_{N-1}|Xi, i~=k)
R2_rest = zeros(N,1);
for i = 1:N
    X_rest = X;
    X_rest(:,i) = [];
    mdl = fitlm(X_rest,y);
    R2_rest(i) = mdl.Rsquared.Ordinary;
end
R2_marg = R2_full - R2_rest;
%calculate the general dominance (Shapley) values
Dominance_value = (R2_marg + R2_indi)/2;
%calculate the standardized domiance (Shapley) values in percentage
Dominance_percentage = Dominance_value./sum(Dominance_value)*100;
%rank the numbers
[~,p] = sort(Dominance_percentage,'descend');
Rank = 1:length(Dominance_percentage);
Rank = Rank';
Rank(p) = Rank;
%Wrap up results in a table
Variable_name = cell(N,1);
for i = 1:N
    Variable_name{i} = strcat('X',num2str(i));
end
tab = table(Variable_name,Dominance_value, Dominance_percentage,Rank); 
end