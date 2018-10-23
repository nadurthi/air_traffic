function index = GenerateIndex(ND,numbasis)
%
% This function computes all permutations of 1-D basis functions.
%
index = (1:numbasis(1))';  %short for canonical_0 - first dimension's nodes: this will be loooped through the dimensions
for ct = 2:ND    %good loop! - over the dimensions
    repel = index; %REPetition-ELement
    repsize = length(index(:,1));  %REPetition SIZE
    repwith = ones(repsize,1);  %REPeat WITH this structure: initialization
    for rs = 2:numbasis(ct)
        repwith = [repwith; ones(repsize,1)*rs];    %update REPeating structure
    end
    index = [repmat(repel,numbasis(ct),1), repwith];    %update canon0
end
index = index;