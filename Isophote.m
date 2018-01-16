%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des isophotes
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

function I = Isophote(N,L,c,S)

I = [];
[n,m,~] = size(N);

disp(['taille de S = ' num2str(size(S))]);

for i = 1:n
    for j = 1:m
        tmp = [N(i,j,1),N(i,j,2),N(i,j,3)];
%         disp(['taille de N(i,j,:) = ' num2str(size(N(i,j,:)))]);
%         disp(['taille de L = ' num2str(size(L))]);
%         disp(['taille de tmp = ' num2str(size(tmp))]);
        I(i,j) = dot(tmp,L);
    end
end
        
res = zeros(n,m,3);

for k = 1:length(c)
    for i = 1:n
        for j = 1:m
            if abs(I(i,j)-c(k)) < 0.01
%             res(i,j) = [S(i,j,1),S(i,j,2),S(i,j,3)];
                res(k,i,j,:) = S(i,j,:);
            end
        end
    end
end

disp(['taille de res = ' num2str(size(res))]);

x = res(:,:,1);
y = res(:,:,2);
z = res(:,:,3);
plot3(x,y,z);



