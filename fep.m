function [K, U]=fep(E,nu,sigma,r, p,t,elemType,gaussType,gaussDeg,print)
%FE 2D function for finding displacements by Alex P. Martinez
%Preprocessor and processor for displacements of hole/s in plate
%Function arguments (all in SI units)
%E:         elastic modulus
%nu:        Poisson's ratio
%Sigma:     applied stress
%p:         node coordinates
%t:         element nodes
%elemType:  element type ('T3', 'T6', 'Q4', 'Q9')
%gaussType: gauss points for standard square or triangle
%(Either 'GAUSS' or 'TRIANGULAR')
%gaussDeg:  quadrature order
%print:     1 to print preprocessing and K info, 0 otherwise
%Outputs:
%K:         global stiffness matrix with boundary conditions
%U:         displacements at nodes


%XXXXXXXXXXXXXXXXXXXXXXX     PREPROCESSOR         XXXXXXXXXXXXXXXXXXXXXXXXX

%BUILD INITIAL MATRICES
%Get element and node numbers and print
nNodes    = size(p,1);
nElements = size(t,1);
nPerElem  = size(t,2);
%Get boundary nodes with certain y values 
%(yBoundary - tol < yNode < yBoundary + tol)
%Boundary nodes
%Tolerance used in boundary selection and y coordinate of boundaries
btol = 1e-4; ybFix = -1; ybForce = 1;
%Fixed nodes
fixedbNodes = find(abs(p(:,2)-ybFix)   < btol);
%Force nodes
forcebNodes = find(abs(p(:,2)-ybForce) < btol);
%Calculate D matrix used in all elements (D matrix for plane stress)
D = E/(1-nu^2)*[1 ,nu,0; nu,1 ,0; 0 ,0 ,(1-nu)/2];
%Define displacement, force vectors and 3d ke matrix
%(MATLAB is faster solving for K*U=F if K is specified as sparce)
U  = zeros(2*nNodes,1); F=U;
K  = sparse(2*nNodes,2*nNodes);
gL = zeros(2*nPerElem,1);


%XXXXXXXXXXXXXXXXXXXXXXX     PROCESSOR            XXXXXXXXXXXXXXXXXXXXXXXXX

%GENERATE K MATRIX
%Get quadrature points in 2D (first number is order of entries in k^e) 
[qW,qP] = quadrature(gaussDeg,gaussType,2);
%Generate ke matrices and put directly into K matrix
tic
for e = 1:nElements
    %k quadrature loop
    for q = 1:size(qW,1) 
        %Calculate derivatives of shape functions
        [N,dNdxi] = lagrange_basis(elemType,qP(q,:));
        elemNodes = p(t(e,:)',:);
        J         = (elemNodes'*dNdxi)';
        dNdx      = dNdxi/J';
        %Make B matrix for element
        B = zeros(3,2*nPerElem);
        B(1,1:2:(2*nPerElem-1)) = dNdx(:,1)';
        B(2,2:2:2*nPerElem    ) = dNdx(:,2)';
        B(3,1:2:(2*nPerElem-1)) = dNdx(:,2)';
        B(3,2:2:2*nPerElem    ) = dNdx(:,1)';
        %Generate K matrix
        %Global location vector
        gL(2*(1:nPerElem)-1) = 2.*t(e,:)-1;
        gL(2*(1:nPerElem))   = 2.*t(e,:);
        K(gL,gL) = K(gL,gL) + B'*D*B*qW(q)*det(J);
    end
end
toc
%GENERATE EXTERNAL FORCE VECTORS
%Sort force nodes based on x value
[fbNxs,indices] = sort(p(forcebNodes,1));
%Calculate lengths between nodes
numfNs = length(forcebNodes);
ls     = diff(fbNxs);
nodalForces = zeros(1,numfNs);
%Choose between linear or quadratic elements
switch elemType
    case {'T3','Q4'}
        %Calculate loads in nodes (consistent loading) adding forces from sides
        %(linear case)
        n = 2:(numfNs-1);
        nodalForces(1)      = ls(1)/2; 
        nodalForces(n)      = (ls(n-1) + ls(n))/2;
        nodalForces(numfNs) = ls(numfNs-1)/2;
    case {'T6','Q8'}
        %Calculate loads in nodes (consistent loading) adding forces from sides
        %(quadratic case)
        n = 3:2:(numfNs-2);
        nodalForces(1)      = qholder(1, 1, ls);
        nodalForces(2)      = qholder(2, 1, ls);
        nodalForces(n)      = qholder(3, n-2, ls) + qholder(1, n, ls);
        nodalForces(n+1)    = qholder(2, n, ls);
        nodalForces(numfNs) = qholder(3, numfNs-2, ls);    
end
%Place forces in force vector
F(forcebNodes(indices)*2) = nodalForces.*sigma;

%APPLY BOUNDARY CONDITIONS
%Zero rows and columns on y component of fixed nodes in K matrix
%Zero leftmost node in x component to prevent free body motion in x axis
%Not needed for F vector since corresponding positions are already zero
[dummy,index] = min(p(fixedbNodes,1));
K(fixedbNodes(index)*2-1,:) = 0; K(:,fixedbNodes(index)*2-1) = 0;
K(fixedbNodes*2,:)          = 0; K(:,fixedbNodes*2)          = 0;
%Put ones in diagonals corresponding to fixed components of nodes
K(fixedbNodes(index)*2-1,fixedbNodes(index)*2-1) = 1;
K(fixedbNodes*2,fixedbNodes*2) = speye(length(fixedbNodes));

%SOLVE SYSTEM
U = K\F; %Same as inv(K)*F in MATLAB notation

%DISPLAY
if print == 1
    %Print processing info
    disp(' ')
    disp('General pre-processing information')
    disp(sprintf(' - Width of square plate (m):    1.0'))
    disp(sprintf(' - Radius of hole (m):           %1.1e',      r))
    disp(sprintf(' - Elastic modulus (Pa):         %1.1e',      E))
    disp(sprintf(' - Poisson''s ratio:              %1.1e',    nu))
    disp(sprintf(' - Applied stress (Pa):          %1.1e',  sigma))
    disp(sprintf(' - Element type:                 %s',  elemType))
    
    disp(sprintf('\n - Number of nodes:              %d',    nNodes))
    disp(sprintf(' - Number of elements:           %d', nElements))
    %Properties of K
    [a,b] = find(K);
    disp(sprintf(' - Non-zero entries in K:        %d',        nnz(K)))
    disp(sprintf(' - Bandwidth of K:               %d', max(abs(a-b))))
    disp(' ')
end