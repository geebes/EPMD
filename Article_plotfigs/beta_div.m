function D2 = beta_div(ZI,ZJ)

% This is effectively presence-absence bray-curtis
D2=zeros(size(ZI'));

a = sum(ZI>0  & ZJ>0 ,2); % n species shared between sites
b = sum(ZI>0  & ZJ==0,2); % n species unique to site 1
c = sum(ZI==0 & ZJ>0 ,2); % n species unique to site 2

D2 = (b + c) ./ (2.*a + b + c);

% With number of cells this is bray-curtis
% D2=zeros(size(ZI'));
% 
% a = sum(min(ZI,ZJ),2); % n cells shared between sites
% b = sum(max(ZI-ZJ,0),2); % n cells unique to site 1
% c = sum(max(ZJ-ZI,0),2); % n cells unique to site 2
% 
% D2 = (b + c) ./ (2.*a + b + c);