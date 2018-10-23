function [xint,wint] = get_colocation(Ninteg, xl, xu)

% Ninteg = [5 5];
% xl = [-10 -1];
% xu = [10 1];


ND = length(xl);

for ct = 1:ND
    [xint, wint] = lgwt(Ninteg(ct),xl(ct),xu(ct)); % Generate Collocation Points
    xinteg(:,ct) = xint(:);
    winteg(:,ct) = wint(:);
end

clear xint wint

[index,xint,wint] = Pcomb(ND,Ninteg,xinteg,winteg);
wint=wint/sum(wint);
xint = xint;