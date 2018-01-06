%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul du curvature plot
%
% Input:
%  - matrix N de num_nxnum_n points de dim 3
%     chaque point a 3 coordonnees (x,y,z)
%     taille de N: num_nxnum_nx3
%  - L dimension des rayons lumineux parallèles
%  - c valeur constante
%  - matrix S des valeurs de l'évaluation du patch de Bézier

% Output:
%  - matrix I avec la grille de |u|x|v| points 3D sur la surface
%    Taille de I: num_nxnum_n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = CurvaturePlot(B,N,u,v)

H = DeuxFormeFond(B,N,u,v);
G = invPremFormeFond(B,u,v);
K = [];

for x = 1:length(u)
    for y = 1:length(v)
        A = mtimes(H(x,y,:),G(x,y,:));
        [~,D] = eig(A);
        K(x,y) = det(D);
    end
end