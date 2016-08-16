function graphSs(elemType, p, t, S)
%This function graphs the following:
%Sparcity graph of K with boundary conditions
%Mesh defined by p and t
%X and y components of displacements at nodes

switch elemType
    case {'T3', 'Q4'}
        %Graph mesh
        %figure();
        %trimesh(t,p(:,1),p(:,2),zeros(size(p,1),1),'EdgeColor',[0,0,0.5])
        %view(2),axis equal,axis off,drawnow

        %Graph x and y components of stress
        figure(); trisurf(t,p(:,1),p(:,2),S)
        axis([-1 1 -1 1 ])
        colormap winter;
        xlabel('x-component of position in plate (m)','FontSize',12)
        ylabel('y-component of position in plate (m)','FontSize',12)
        zlabel('y-component stress (Pa)','FontSize',12)
        title('Y-component of stress as a function of position (m)','FontSize',14)
    
    case 'T6'
        %Graph mesh
        %figure();
        %trimesh(t(:,[1 2 3]),p(:,1),p(:,2),zeros(size(p,1),1),'EdgeColor',[0,0,0.5])
        %view(2),axis equal,axis off,drawnow
        
        %Graph x and y components of displacement
        figure(); trisurf(t(:,[1 2 3]),p(:,1),p(:,2),S)
        axis([-1 1 -1 1 ])
        colormap winter;
        
    case 'Q8'
        %Graph mesh
        %figure();
        %trimesh(t(:,[1 2 3 4]),p(:,1),p(:,2),zeros(size(p,1),1),'EdgeColor',[0,0,0.5])
        %view(2),axis equal,axis off,drawnow
        
        %Graph x and y components of displacement
        figure(); trisurf(t(:,[1 2 3 4]),p(:,1),p(:,2),S)
        axis([-1 1 -1 1 ])
        colormap winter;
end