function graphDs(elemType, p, t, K, U)
%This function graphs the following:
%Sparcity graph of K with boundary conditions
%Mesh defined by p and t
%X and y components of displacements at nodes


%Show sparcity graph in new window, in black ('k') and with marker size 5
%(with zeroed columns and rows)
figure();
spy(K,'k',4)
title('Sparcity graph of stiffness matrix K','FontSize',14)

switch elemType
    case {'T3', 'Q4'}
        %Graph mesh
        %figure();
        %trimesh(t,p(:,1),p(:,2),zeros(size(p,1),1),'EdgeColor',[0,0,0.5])
        %view(2),axis equal,axis off,drawnow

        %Graph x and y components of displacement
        n = 2:2:size(p,1)*2;
        figure(); trisurf(t,p(:,1),p(:,2),U(n))
        axis([-1 1 -1 1 ])
        colormap winter;
        xlabel('x-component of position in plate (m)','FontSize',12)
        ylabel('y-component of position in plate (m)','FontSize',12)
        zlabel('y-component of node dispalcements (m)','FontSize',12)
        title('Y-component of node displacements as a function of position (m)','FontSize',14)
%         figure(); trisurf(t,p(:,1),p(:,2),U(n-1))
%         axis([-1 1 -1 1 ])
%         colormap winter;
    
    case 'T6'
        %Graph mesh
        %figure();
        %trimesh(t(:,[1 2 3]),p(:,1),p(:,2),zeros(size(p,1),1),'EdgeColor',[0,0,0.5])
        %view(2),axis equal,axis off,drawnow
        
        %Graph x and y components of displacement
        n = 2:2:size(p,1)*2;
        figure(); trisurf(t(:,[1 2 3]),p(:,1),p(:,2),U(n))
        axis([-1 1 -1 1 ])
        colormap winter;
        %figure(); trisurf(t(:,[1 2 3]),p(:,1),p(:,2),U(n-1))
        %axis([-1 1 -1 1 ])
        %colormap winter;
        
    case 'Q8'
        %Graph mesh
        %figure();
        %trimesh(t(:,[1 2 3 4]),p(:,1),p(:,2),zeros(size(p,1),1),'EdgeColor',[0,0,0.5])
        %view(2),axis equal,axis off,drawnow
        
        %Graph x and y components of displacement
        n = 2:2:size(p,1)*2;
        figure(); trisurf(t(:,[1 2 3 4]),p(:,1),p(:,2),U(n))
        axis([-1 1 -1 1 ])
        colormap winter;
        figure(); trisurf(t(:,[1 2 3 4]),p(:,1),p(:,2),U(n-1))
        axis([-1 1 -1 1 ])
        colormap winter;
end