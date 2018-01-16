%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul du curvature plot avec la courbure de Gauss
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

function K = CurvaturePlot(B,N,S,u,v,np)

H = DeuxFormeFond(B,N,u,v);
G = InvPremFormeFond(B,u,v);
K = [];

for x = 1:length(u)
    for y = 1:length(v)
        K(x,y) = (((H(x,y)/G(x,y)) + 1)/2)*255;
    end
end