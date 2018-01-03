%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des isophotes
%
% Input:
%  - matrix N de num_nxnum_n points de dim 3
%     chaque point a 3 coordonnees (x,y,z)
%     taille de N: num_nxnum_nx3
%  - L dimension des rayons lumineux parallèles
%  - c valeur constante

% Output:
%  - matrix I avec la grille de |u|x|v| points 3D sur la surface
%    Taille de I: num_nxnum_n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = Isophote(N,L,c,S)

I = [];

for i = 1:size(N,1)
    for j = 1:size(N,2)
        I(i,j) = dot(N(i,j,:),L);
    end
end
        
res = [];

for i = 1:size(N,1)
    for j = 1:size(N,2)
        if norm(I(i,j)-c) < 0.001
            res(i,j) = S(i,j,:);
        else 
            res(i,j) = [0,0,0];
        end
    end
end

plot3(res);



