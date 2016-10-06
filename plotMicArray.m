function plotMicArray(mic_dirs_deg, R)
% PLOTMICARRAY plots the arrangement of the microphones in a spherical array
%
%   mic_dirs:   the directions of the microphones in [azi1 elev1; ...]
%       convention
%   R:          radius of the array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PLOTMICARRAY.M - 11/7/2013
% Archontis Politis, archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mic_dirs_rad = mic_dirs_deg*pi/180;
Nmic = size(mic_dirs_deg,1);

hold on
% set up unit sphere information
numSphereFaces = 20;
[unitSphereX, unitSphereY, unitSphereZ] = sphere(numSphereFaces);
% radius of each sphere
spheresRadius = ones(Nmic,1)*0.05*R;
%plot 3d axes
line([0 2*R],[0 0], [0 0],'color','r');
text(2*R,0,0,'x','Color','r','FontSize',24);
line([0 0],[0 2*R], [0 0],'color','g');
text(0,2*R,0,'y','Color','g','FontSize',24);
line([0 0],[0 0], [0 2*R],'color','b');
text(0,0,2*R,'z','Color','b','FontSize',24);

spheresX = R*cos(mic_dirs_rad(:,1)).*cos(mic_dirs_rad(:,2));
spheresY = R*sin(mic_dirs_rad(:,1)).*cos(mic_dirs_rad(:,2));
spheresZ = R*sin(mic_dirs_rad(:,2));
% for each given sphere, shift the scaled unit sphere by the
% location of the sphere and plot
for i=1:Nmic
    sphereX = spheresX(i) + unitSphereX*spheresRadius(i);
    sphereY = spheresY(i) + unitSphereY*spheresRadius(i);
    sphereZ = spheresZ(i) + unitSphereZ*spheresRadius(i);
    h = surf(sphereX, sphereY, sphereZ);
    set(h, 'FaceColor','b', 'FaceLighting', 'gouraud')
    
    text(1.1*spheresX(i),1.1*spheresY(i),1.1*spheresZ(i), num2str(i),'color','w','FontSize',24)    
end
sphereX = unitSphereX*R;
sphereY = unitSphereY*R;
sphereZ = unitSphereZ*R;
h = surf(sphereX, sphereY, sphereZ);
set(h, 'FaceColor','c', 'FaceLighting', 'gouraud')
light('Position',[0 0 1],'Style','infinite');
material shiny
hold off
axis equal
grid on
view(-30,15)
set(gca,'visible','off')
