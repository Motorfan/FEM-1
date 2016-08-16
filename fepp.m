function [Sx, Sy, Sxy] = fepp(E,nu,p,t,U,elemType,gaussType,print)
%FE 2D function for finding stresses from displacements by Alex P. Martinez
%Postprocessor for displacements of hole/s in plate
%Function arguments (all in SI units)
%E:         elastic modulus
%nu:        Poisson's ratio
%p:         node coordinates
%t:         element nodes
%elemType:  element type ('T3', 'T6', 'Q4', 'Q9')
%gaussType: gauss points for standard square or triangle
%(Either 'GAUSS' or 'TRIANGULAR')
%gaussDeg:  quadrature order
%print:     1 to print postprocessing info, 0 otherwise
%Outputs:
%S:         stresses at nodes
%diffS:     difference between lowest and highest stresses 
%in nodes calculated at each element before averaging 


%XXXXXXXXXXXXXXXXXXXXXXX     POSTPROCESSOR         XXXXXXXXXXXXXXXXXXXXXXXX

%Calculate D matrix used in all elements (D matrix for plane stress)
D = E/(1-nu^2)*[1 ,nu,0; nu,1 ,0; 0 ,0 ,(1-nu)/2];
%Get element and node numbers
nNodes    = size(p,1);
nElements = size(t,1);
nPerElem  = size(t,2);


%CALCULATE STRESSES AT GAUSS POINTS
%Get quadrature points in 2D (first number is order of entries in k^e) 
switch elemType
    case 'T3' 
        gaussDeg = 1;
        order = 1;
    case 'T6'
        gaussDeg = 2;
        order = (1:3)';
    case {'Q4','Q8'}
        gaussDeg = 2;
        order = [4;2;1;3];
end     
[qW,qP] = quadrature(gaussDeg,gaussType,2);
sizeQuad  = size(qW,1);
orderSize = size(order,1);

%Declare stress at quad points arrays
Sqx  = zeros(orderSize, orderSize);
Sqy  = zeros(orderSize, orderSize);
Sqxy = zeros(orderSize, orderSize);
adj  = zeros(nNodes, 1);
Sx   = zeros(nNodes, 1);
Sy   = zeros(nNodes, 1);
Sxy  = zeros(nNodes, 1);
T    = zeros(nPerElem, orderSize);

%Stress calculation
for e = 1:nElements
    %k quadrature loop
    for q = 1:orderSize
        %Calculate derivatives of shape functions
        [N,dNdxi] = lagrange_basis(elemType,qP(order(q),:));
        elemNodes = p(t(e,:)',:);
        J         = (elemNodes'*dNdxi)';
        invJt     = inv(J');
        dNdx      = dNdxi*invJt;
        %Make B matrix for element
        B = zeros(3,2*nPerElem);
        B(1,1:2:(2*nPerElem-1)) = dNdx(:,1)';
        B(2,2:2:2*nPerElem    ) = dNdx(:,2)';
        B(3,1:2:(2*nPerElem-1)) = dNdx(:,2)';
        B(3,2:2:2*nPerElem    ) = dNdx(:,1)';
        %Actual stress calculation
        %Global location vector
        gL(2*(1:nPerElem)-1) = 2.*t(e,:)-1;
        gL(2*(1:nPerElem))   = 2.*t(e,:);
        holder = D*B*U(gL);
        Sqx(e,q) = holder(1); Sqy(e,q) = holder(2);
        Sqxy(e,q) = holder(3);
    end
end


%EXTRAPOLATE POINTS
%Decide how to extrapolate based on element type
switch elemType
    case 'T3' 
        for e = 1:nElements
            adj(t(e,:),1) = adj(t(e,:),1) + 1;
            
            Sx(t(e,:), 1) = Sx(t(e,:), 1) + Sqx(e,1);
            Sy(t(e,:), 1) = Sy(t(e,:), 1) + Sqy(e,1);
            Sxy(t(e,:),1) = Sxy(t(e,:),1) + Sqxy(e,1);
        end

    case 'T6'
        elemType = 'T3';
        Ns = [0,0; 1,0; 0,1;  0.5,0; 0.5,0.5; 0,0.5];
        
        for q = 1:nPerElem
            newP = (Ns(q,:) - qP(1,1))./(qP(2,1)-qP(1,1)); 
            
            [N,dNdxi] = lagrange_basis(elemType, newP);
            T(q,:) = N';
        end
        
        for e = 1:nElements
            adj(t(e,:),1) = adj(t(e,:),1) + 1;
            
            Sx(t(e,:),1)   = Sx(t(e,:),1)   + T*Sqx(e,:)';
            Sy(t(e,:),1)   = Sy(t(e,:),1)   + T*Sqy(e,:)';
            Sxy(t(e,:),1)  = Sxy(t(e,:),1)  + T*Sqxy(e,:)';
        end
        
    case 'Q4'
        for q = 1:orderSize
            [N,dNdxi] = lagrange_basis(elemType,1./qP(order(q),:));
            T(q,:) = N';
        end
        
        for e = 1:nElements
            adj(t(e,:),1) = adj(t(e,:),1) + 1;
            
            Sx(t(e,:), 1) = Sx(t(e,:), 1) + T*Sqx(e, :)';
            Sy(t(e,:), 1) = Sy(t(e,:), 1) + T*Sqy(e, :)';
            Sxy(t(e,:),1) = Sxy(t(e,:),1) + T*Sqxy(e,:)';
        end
        
    case 'Q8'
        elemType = 'Q4';
        Ns = [-1,-1; 1,-1; 1,1; -1,1;  0,-1; 1,0; 0,1; -1,0];
        
        for q = 1:nPerElem
            [N,dNdxi] = lagrange_basis(elemType, Ns(q,:)./qP(1,1));
            T(q,:) = N';
        end
        
        for e = 1:nElements
            adj(t(e,:),1) = adj(t(e,:),1) + 1;
            
            Sx(t(e,:),1)   = Sx(t(e,:),1)   + T*Sqx(e,:)';
            Sy(t(e,:),1)   = Sy(t(e,:),1)   + T*Sqy(e,:)';
            Sxy(t(e,:),1)  = Sxy(t(e,:),1)  + T*Sqxy(e,:)';
        end
end 
Sx  = Sx./adj;
Sy  = Sy./adj;
Sxy = Sxy./adj;


%DISPLAY
if print == 1
    %Print postprocessing info
    [maxSx, I]  = max(abs(Sx));
    maxSx  = sign(Sx(I))*maxSx;
    [maxSy, I]  = max(abs(Sy));
    maxSy  = sign(Sy(I))*maxSy;
    [maxSxy, I] = max(abs(Sxy));
    maxSxy = sign(Sxy(I))*maxSxy;
    disp('General post-processing information')
    disp(sprintf(' - Max mag x  stress component:  %8.2f', maxSx))
    disp(sprintf(' - Max mag y  stress component:  %8.2f', maxSy))
    disp(sprintf(' - Max mag xy stress component:  %8.2f', maxSxy))
    disp(' ')
end