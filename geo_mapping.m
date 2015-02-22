function [x_plus] = geo_mapping(filename,d_alpha1,d_alpha2)
%  Computes a geometric mapping for the two ellipse flow
%  problem given \delta\alpha_1 and \delta\alpha_2.
% 
%  The mapping is created by solving the two-dimensional
%  elasticity problem with prescribed parametric perturbations
%  as Dirichlet boundary conditions.
%%

  addpath('../fem_functions')  % keep it relative

  % Elastic constants
  material.young = 1e5;  
  material.nu    = 0.3;

  % Solution parameter
  Penalty=10^20;

  %% Load mesh
  [x,e_conn,boundary] = read_msh(filename);
    
  % Renumber unknowns to reduce bandwidth
  stream = RandStream.getGlobalStream;
  reset(stream);
  
  [x,e_conn,p_inv] = tri_mesh_rcm(x,e_conn);
  for n=1:length(boundary)
    boundary(n).nodes = p_inv(boundary(n).nodes);
  end
  
  if ( strcmp(boundary(5).name,'all') )
    ibcd = boundary(5).nodes';
  else
    error('problem with boundary labels')
  end

  n_nodes=size(x,1); n_dof=2*n_nodes;

  %% Initialize Dirichlet boundary conditions to zero
  %  Calculate global indices of Dirichlet conditions
  n_dbc        = 2*length(ibcd);
  ibc          = zeros(n_dbc,1); 
  ibc(1:2:end) = 2*ibcd-1; 
  ibc(2:2:end) = 2*ibcd;

  ubc = zeros(n_dof,1);     % could be sparse if too large


  % slightly rotate the upper ellipse
%   a = 1;
%   b = 0.5;
%   alpha1 = pi/6;
  xc = 0;
  yc = 2;
  rot = [ cos(d_alpha1) -sin(d_alpha1); sin(d_alpha1) cos(d_alpha1) ];

  for n=boundary(3).nodes
    c = [ x(n,1)-xc; x(n,2)-yc ];
    d = rot*c;
   
    ubc(2*n-1) = d(1)-c(1);
    ubc(2*n  ) = d(2)-c(2);
  end
 

  % slightly rotate the lower ellipse
%   a = 1; %#ok
%   b = 0.5;
%   alpha2 = pi/4;

  xc = 2;
  yc = -1;
  rot = [ cos(d_alpha2) -sin(d_alpha2); sin(d_alpha2) cos(d_alpha2) ];

  for n=boundary(4).nodes
    c = [ x(n,1)-xc; x(n,2)-yc ];
    d = rot*c;
   
    ubc(2*n-1) = d(1)-c(1);
    ubc(2*n  ) = d(2)-c(2);
  end


  %% Integrate the stiffness matrix and incorporate body forces 
  [K,f] = elasticity_matrices(x,e_conn,material);

  % Apply boundary conditions by penalizing the stiffness
  % matrix and augmenting the right-hand side
  K(ibc,ibc) = K(ibc,ibc) + Penalty*speye(n_dbc);
  f          =          f + Penalty*ubc;

  % Solve to compute all displacements
  u=K\f;

  u = reshape(u',2,n_nodes)';
  x_plus = x + u;

  
  %% Optional post-processing
  % Compute the shear energy density & von Mises effective stress
%   
%   twod_to_vtk('elastic.vtk',x_plus,e_conn,[shear_energy vonMises],u,...
%                               {'shearEnergy','vonMises','displacement'})
                   
end
