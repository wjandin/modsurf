%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul de la deuxieme forme fondamentale
%
% Input:
%  - matrix B de 16 points de controle de dim 3
%     chaque point de controle a 3 coordonnees (x,y,z)
%     taille de B: 4x4x3
%     B(:,:,k) keme coordonnee des 16  points de controle, k=1,2,3
%     B(i,j,:) les 3 coordonnees du point de controle b_ij
%  - matrix N du champ des normales
%  - u Vecteur de |u|=length(u) valeurs de parametre en u
%  - v Vecteur de |v|=length(v) valeurs de parametre en v

% Output:
%  - matrix H avec la grille de |u|x|v| points 3D sur la surface
%    La structure de H est similaire a celle de B.
%    Taille de H: |u|x|v|x2x2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H = DeuxFormeFond(B,N,u,v)

H = [];

for x = 1:length(u)
    for y = 1:length(v)
        
        X_uu = [0 0 0];
        for i = 1:2
            for j = 1:4
                X_uu(1) = X_uu(1) + (B(i+2,j,1)-2*B(i+1,j,1)+B(i,j,1)) * nchoosek(1,i-1) * power(1-u(x),1-i+1) * power(u(x),i-1) * nchoosek(3,j-1) * power(1-v(y),3-j+1) * power(v(y),j-1);
                X_uu(2) = X_uu(2) + (B(i+2,j,2)-2*B(i+1,j,2)+B(i,j,2)) * nchoosek(1,i-1) * power(1-u(x),1-i+1) * power(u(x),i-1) * nchoosek(3,j-1) * power(1-v(y),3-j+1) * power(v(y),j-1);
                X_uu(3) = X_uu(3) + (B(i+2,j,3)-2*B(i+1,j,3)+B(i,j,3)) * nchoosek(1,i-1) * power(1-u(x),1-i+1) * power(u(x),i-1) * nchoosek(3,j-1) * power(1-v(y),3-j+1) * power(v(y),j-1);
            end    
        end
        X_uu = X_uu * 6;
        
        
        X_vv = [0 0 0];
        for i = 1:4
            for j = 1:2
                X_vv(1) = X_vv(1) + (B(i,j+2,1)-2*B(i,j+1,1)+B(i,j,1)) * nchoosek(3,i-1) * power(1-u(x),3-i+1) * power(u(x),i-1) * nchoosek(1,j-1) * power(1-v(y),1-j+1) * power(v(y),j-1);
                X_vv(2) = X_vv(2) + (B(i,j+2,2)-2*B(i,j+1,2)+B(i,j,2)) * nchoosek(3,i-1) * power(1-u(x),3-i+1) * power(u(x),i-1) * nchoosek(1,j-1) * power(1-v(y),1-j+1) * power(v(y),j-1);
                X_vv(3) = X_vv(3) + (B(i,j+2,3)-2*B(i,j+1,3)+B(i,j,3)) * nchoosek(3,i-1) * power(1-u(x),3-i+1) * power(u(x),i-1) * nchoosek(1,j-1) * power(1-v(y),1-j+1) * power(v(y),j-1);
            end    
        end
        X_vv = X_vv * 6;
        
        X_uv = [0 0 0];
        for i = 1:3
            for j = 1:3
                X_uv(1) = X_uv(1) + (B(i+1,j+1,1)-B(i,j+1,1)-B(i+1,j,1)+B(i,j,1)) * nchoosek(2,i-1) * power(1-u(x),2-i+1) * power(u(x),i-1) * nchoosek(2,j-1) * power(1-v(y),2-j+1) * power(v(y),j-1);
                X_uv(2) = X_uv(2) + (B(i+1,j+1,2)-B(i,j+1,2)-B(i+1,j,2)+B(i,j,2)) * nchoosek(2,i-1) * power(1-u(x),2-i+1) * power(u(x),i-1) * nchoosek(2,j-1) * power(1-v(y),2-j+1) * power(v(y),j-1);
                X_uv(3) = X_uv(3) + (B(i+1,j+1,3)-B(i,j+1,3)-B(i+1,j,3)+B(i,j,3)) * nchoosek(2,i-1) * power(1-u(x),2-i+1) * power(u(x),i-1) * nchoosek(2,j-1) * power(1-v(y),2-j+1) * power(v(y),j-1);
            end    
        end
        X_uv = X_uv * 9;

        %disp(['taille de N : ' num2str(size(N(x,y,:)))]);
        tmp = [N(i,j,1),N(i,j,2),N(i,j,3)];
        H12 = dot(X_uv,tmp);
        H11 = dot(X_uu,tmp);
        H22 = dot(X_vv,tmp);
        H(x,y,:) = H11 * H22 - H12 * H12;
    end
end