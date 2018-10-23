function [index,xint,wint] = Pcomb(ND,numquad,xinteg,winteg)
index = (1:numquad(1))';  %short for canonical_0 - first dimension's nodes: this will be loooped through the dimensions
for ct = 2:ND    %good loop! - over the dimensions
    repel = index; %REPetition-ELement
    repsize = length(index(:,1));  %REPetition SIZE
    repwith = ones(repsize,1);  %REPeat WITH this structure: initialization
    for rs = 2:numquad(ct)
        repwith = [repwith; ones(repsize,1)*rs];    %update REPeating structure
    end
    index = [repmat(repel,numquad(ct),1), repwith];    %update canon0
end
index = index;

xint = xinteg(index(:,1),1);
wint = winteg(index(:,1),1);
for j = 2:ND %good loop! - run through the number of dimensions!
    xint = [xint, xinteg(index(:,j),j)];
    wint = wint.*winteg(index(:,j),j);
end