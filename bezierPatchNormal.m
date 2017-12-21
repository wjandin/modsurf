%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul du champs de normales aux parametres u,v
%
% Input:
%  - matrix B de 16 points de controle de dim 3
%     chaque point de controle a 3 coordonnees (x,y,z)
%     taille de B: 4x4x3
%     B(:,:,k) keme coordonnee des 16  points de controle, k=1,2,3
%     B(i,j,:) les 3 coordonnees du point de controle b_ij
%  - u Vecteur de |u|=length(u) valeurs de parametre en u
%  - v Vecteur de |v|=length(v) valeurs de parametre en v

% Output:
%  - matrix S avec la grille de |u|x|v| points 3D sur la surface
%    La structure de S est similaire a celle de B.
%    Taille de S: |u|x|v|x3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = bezierPatchNormal(B,u,v)

N=[];

for x = 1:length(u)
    for y = 1:length(v)
        
        X_u = [0 0 0];
        for i = 1:3
            for j = 1:4
                X_u(1) = X_u(1) + (B(i+1,j,1)-B(i,j,1)) * nchoosek(2,i-1) * power(1-u(x),2-i+1) * power(u(x),i-1) * nchoosek(3,j-1) * power(1-v(y),3-j+1) * power(v(y),j-1);
                X_u(2) = X_u(2) + (B(i+1,j,2)-B(i,j,2)) * nchoosek(2,i-1) * power(1-u(x),2-i+1) * power(u(x),i-1) * nchoosek(3,j-1) * power(1-v(y),3-j+1) * power(v(y),j-1);
                X_u(3) = X_u(3) + (B(i+1,j,3)-B(i,j,3)) * nchoosek(2,i-1) * power(1-u(x),2-i+1) * power(u(x),i-1) * nchoosek(3,j-1) * power(1-v(y),3-j+1) * power(v(y),j-1);
            end    
        end
        X_u = X_u * 3;
        
        
        X_v = [0 0 0];
        for i = 1:4
            for j = 1:3
                X_v(1) = X_v(1) + (B(i,j+1,1)-B(i,j,1)) * nchoosek(3,i-1) * power(1-u(x),3-i+1) * power(u(x),i-1) * nchoosek(2,j-1) * power(1-v(y),2-j+1) * power(v(y),j-1);
                X_v(2) = X_v(2) + (B(i,j+1,2)-B(i,j,2)) * nchoosek(3,i-1) * power(1-u(x),3-i+1) * power(u(x),i-1) * nchoosek(2,j-1) * power(1-v(y),2-j+1) * power(v(y),j-1);
                X_v(3) = X_v(3) + (B(i,j+1,3)-B(i,j,3)) * nchoosek(3,i-1) * power(1-u(x),3-i+1) * power(u(x),i-1) * nchoosek(2,j-1) * power(1-v(y),2-j+1) * power(v(y),j-1);
            end    
        end
        X_v = X_v * 3;

        N(x,y,:) = cross(X_u,X_v) / norm(cross(X_u,X_v));
    end
end


