%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des isophotes
%
% Input:
%  - matrix N de num_nxnum_n points de dim 3
%     chaque point a 3 coordonnees (x,y,z)
%     taille de N: num_nxnum_nx3
%  - L dimension des rayons lumineux paralleles
%  - c valeur constante
%  - matrix S des valeurs de l'evaluation du patch de Bezier

% Output:
%  - matrix I avec la grille de |u|x|v| points 3D sur la surface
%    Taille de I: num_nxnum_n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = Isophote(N,L)

I = [];
[n,m,~,~] = size(N);
for i = 1:n
    for j = 1:m
        tmp = [N(i,j,1),N(i,j,2),N(i,j,3)];
        I(i,j) = dot(tmp(:),L(:));
    end
end



