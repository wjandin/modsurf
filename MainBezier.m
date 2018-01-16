close all, clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul et visualisation de surfaces de Bezier bi-cubiques
%
% Input: fichier ascii avec 16 x np points de controle
%
% Les patchs sont individuellement evalues en un 
% nombre fixe de parametres (u,v).
%
% Parametres: 
% num_p : nombre de valeurs de parametres = nombre de points de controle
% num_n : nombre de normales calcules (pour le calcul des isophotes)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BezierSurf = load('surface1');  % read control points
%BezierSurf = load('surface2'); % read control points
%BezierSurf = load('surface3'); % read control points
BezierSurf = load('surface4'); % read control points
%load('teapot'); %loading matrix B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_p=100;
%num_p=20;                    % nombre de valeurs de parametre en direction u et v
			     % plus num_p est grand plus la surface paraitra lisse
			     % et plus le calcul sera long
num_n=100;		     % nombre de normales en direction u et en v calcules.

%-------------------------------------------------
[nb,~] = size(BezierSurf) % nombre de points de controle dans le fichier
np = floor(nb/16) % nombre de patches composant la surface
                   % Il faudrait verifier que nb est un multiple de 16
% Matrice B des points de controle
for k=1:np
  for i=1:4
    for j=1:4
      B(i,j,1,k) = BezierSurf((i-1)*4+j+(k-1)*16,1);
      B(i,j,2,k) = BezierSurf((i-1)*4+j+(k-1)*16,2);
      B(i,j,3,k) = BezierSurf((i-1)*4+j+(k-1)*16,3);
    end
  end
end
B

% La matrice B stocke tous les points de controle de tous les patchs
% B(:,:,:,k) sont tous les points de controle du patch k
% La dimension de B(:,:,:,k) est 4 x 4 x 3, i.e., 16 points de controle
% à 3 coordonnees (x,y,z)

% B(:,:,1,k): x-coordonnees des points de controle du patch k comme matrice 4 x 4
% B(:,:,2,k): y-coordonnees des points de controle du patch k comme matrice 4 x 4
% B(:,:,3,k): z-coordonnees des points de controle du patch k comme matrice 4 x 4

% ------------------------------------
% num_p+1 valeurs de parametres uniformes: 0,1,2,...,num_p en u et v
u = linspace(0,1,num_p); 
v = u; 


%  ------------------------------------
% Cubic Bezier patches 
for k=1:np
    S(:,:,:,k)=bezierPatchEval(B(:,:,:,k),u,v); %evaluation du patch k en (num_p+1) x (num_p+1) points
end


% ------------------------------------
% Normal vectors of Cubic Bezier patches 
u=linspace(0,1,num_n); v=u;  %parametrisation uniforme (num_n+1)x (num_n+1) valeurs de parametre
for k=1:np
    N(:,:,:,k)=bezierPatchNormal(B(:,:,:,k),u,v); %vecteurs normaux du patch k
end


% ------------------------------------
% Computing Isophotes
L = [0,0,1];
for k=1:np
    I(:,:,k)=Isophote(N(:,:,:,k),L);
end

c = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];

list_x = [];
list_y = [];
list_z = [];
[n,m,~,~] = size(N);

for z = 1:length(c)
    for k = 1:np
        for i = 1:n
            for j = 1:m
                if abs(I(i,j,k)-c(z)) < 0.01
                    list_x(end+1) = S(i,j,1,k);
                    list_y(end+1) = S(i,j,2,k);
                    list_z(end+1) = S(i,j,3,k);
                end
            end
        end
    end
end

figure, hold on
for k=1:np
    surface(S(:,:,1,k),S(:,:,2,k),S(:,:,3,k))
end
shading interp
title('\bf Surface de Bezier avec Interpolated Shading et isophotes');
view(3); box;  view(21,19)
plot3(list_x,list_y,list_z,'.k');

% ------------------------------------
% Computing Curvature Plot
% for k=1:np
%     K(:,:,k)=CurvaturePlot(B(:,:,:,k),N(:,:,:,k),S(:,:,:,k),u,v,np);
% end
% figure, hold on
% for k=1:np
%     surface(S(:,:,1,k),S(:,:,2,k),S(:,:,3,k),K(:,:,k))
% end
% title('\bf Surface de Bezier avec Curvature Plot');
% view(3); box;  view(21,19)

% ------------------------------------
% Visualisation d'un patch/surface de Bezier
  %plotBezierPatch3D(B(:,:,:,k),S(:,:,:,k)) % plot d'un seul patch k
  %plotBezierSurface3D(B,S)		   % plot de tous les np patches
