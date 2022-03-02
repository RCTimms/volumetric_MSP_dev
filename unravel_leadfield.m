function L_u=unravel_leadfield(L);
% Get first lead field to learn dimensions
naught=L{1};

% Pre-allocate
L_u=zeros(size(naught,1),size(L,1)*3);

% Unpack leadfield
for i=1:size(L,1);
    L_u(:, (3*i - 2):(3*i))  = L{i};
end
end