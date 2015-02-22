function []=peridynamic_1d
%-------------------------------------------------------------------------------
%  peridynamic_1d.m - A simple implementation of a 1D peridynamic model
%
%  Copyright (c) 2014, Jeff Borggaard, Virginia Tech
%  Version: 1.0
%
%  Usage:    [] = peridynamic_1d
%
%  Variables:  
%              elements_mesh
%                        Density of elements in (x_l,x_r)
%
%  Functions:  f_function
%                        interelement force function
%              b_function
%                        external forces
%-------------------------------------------------------------------------------
  addpath('../fem_functions')
  
  %-----------------------------------------------------------------------------
  %  Define "Input" Parameters
  %-----------------------------------------------------------------------------
  elements_mesh    = 80;
  x_l              = 0;
  x_r              = 1;

  ell              = 0.1;     % defines the force horizon
  
  
  element_type     = 'continuous_linear';

  %-----------------------------------------------------------------------------
  %  Setup Geometry
  %-----------------------------------------------------------------------------
  if ( strcmp(element_type,'continuous_linear') )
    xb(:,1)      = [x_l; x_r];
    e_connb(1,:) = [1 2];
    rho(1)       = elements_mesh;
    [x,e_conn]   = oned_mesh(xb,e_connb,rho);
  end

  if ( strcmp(element_type,'continuous_quadratic') )
    xb(:,1)      = [x_l; .5*(x_l+x_r); x_r];
    e_connb(1,:) = [1 2 3];
    rho(1)       = elements_mesh;
    [x,e_conn]   = oned_mesh(xb,e_connb,rho);
  end

  if ( strcmp(element_type,'continuous_cubic') )
    xb(:,1)      = [x_l; (x_l+x_r)/3; 2*(x_l+x_r)/3; x_r];
    e_connb(1,:) = [1 2 3 4];
    rho(1)       = elements_mesh;
    [x,e_conn]   = oned_mesh(xb,e_connb,rho);
  end

  n_nodes = size(x,1);
  ide     = 1:n_nodes;
  
  %-----------------------------------------------------------------------------
  %  Newton Solver
  %-----------------------------------------------------------------------------
  step_tolerance = 1e-6;
  resid_tolerance = 1e-6;

  %-----------------------------------------------------------------------------
  % Provide the initial guess for the displacement, u
  %-----------------------------------------------------------------------------
  u = 2*linspace(0,1,n_nodes)-1;
    
  %  Compute the initial residual
  [R] = weak_resid(x,e_conn,ide,u,ell);

  %  Begin Newton iterations
  converged = false;
  diverged  = false;
  iteration = 0;
  while( ~converged && ~diverged )
    plot(x,u)
    pause(.01)
    
    %  Compute the Jacobian via finite differences (implement sparse fd later)
    step = 1e-5 * norm(u,inf);
    J = zeros(n_nodes,n_nodes);
    for i=1:n_nodes
      up = u;
      up(i) = up(i) + step;
      [Rp] = weak_resid(x,e_conn,ide,up,ell);
      J(:,i) = (Rp-R)/step;
    end
    
%     mid = round((n_nodes+1)/2);
%     J(mid,:) = zeros(1,n_nodes);  J(mid,mid) = 1;
%     R(mid)   = 0;

    J(n_nodes+1,1:n_nodes) = ones(1,n_nodes);
    R(n_nodes+1)           = 0;  %imposes an "average is zero" condition
    
    s = -J\R;
    
    
    iteration = iteration + 1;
    u = u + s';
    
    %  Compute the residual using the current iterate
    [R] = weak_resid(x,e_conn,ide,u,ell);
    
    fprintf('Iteration %d:\n',iteration)
    fprintf('  || step || = %g\n',norm(s))
    fprintf('  || resl || = %g\n\n',norm(R))
    converged = ( norm(s)<step_tolerance ) || ( norm(R)<resid_tolerance );
    diverged  = ( iteration>3 ) || ( norm(R)>1e8 );

  end
  
  figure
  plot(x,u+x','bo')
  
end

%-------------------------------------------------------------------------------
%  Supporting Functions
%-------------------------------------------------------------------------------
function [R] = weak_resid(x,e_conn,ide,u,ell)

  n_gauss     = 4;
  [r,w]       = oned_gauss(n_gauss);
  
  [n_elements,nel_dof] = size(e_conn);
  n_equations          = max(max(ide));

  ide = 1:n_equations;

  %-----------------------------------------------------------------------------
  %  Build the weak residual
  %-----------------------------------------------------------------------------
  R = zeros (n_equations,1);
  
  for n_el=1:n_elements
    % compute value of each test function and their spatial derivaties
    % at the integration points
    nodes_local     = e_conn(n_el,:);
    x_local         = x(nodes_local,:);
    [x_g,w_g,phi,~] = oned_shape(x_local,r,w);

    % compute the value of functions at the Gauss points
    u_local = u(nodes_local);
    
    u_g   = phi*u_local';
    b_g   = b_function(x_g);

    %---------------------------------------------------------------------------
    %  Integrate the weak form of the equations (element contributions)
    %---------------------------------------------------------------------------

    %  Integrate the nonlocal terms
    
    %  Calculate elements in the Horizon of x
    horizon = 1:n_elements;    % for now
    
    int_f = zeros(size(x_g));
    for h_el=horizon
      nodes_local1                = e_conn(h_el,:);
      x_local1                    = x(nodes_local1,:);
      [xi_g,wprime_g,phi_prime,~] = oned_shape(x_local1,r,w);
      
      u_local1 = u(nodes_local1);
      uprime_g = phi_prime*u_local1';   % u at the xi Gauss points
      
      for i=1:n_gauss
        chi_g    = horizon_function(x_g(i),ell,xi_g);
        f_g      = f_function(uprime_g-u_g(i),xi_g-x_g(i));
        int_f(i) = int_f(i) + sum(chi_g .* f_g .* wprime_g);
      end
      
    end
    
    f_loc = oned_f_int( int_f, phi, w_g ) + oned_f_int(b_g, phi, w_g);
    %---------------------------------------------------------------------------
    % Assemble contributions into the global system matrix
    %---------------------------------------------------------------------------
    for n_t=1:nel_dof
      n_test = ide(nodes_local(n_t));
      if (n_test > 0)  % this is an unknown, fill the row
        R(n_test) = R(n_test) + f_loc(n_t);
      end
    end
    
  end
end


function b = b_function(x)
  b = 100*((x>.625&x<.675).*exp(-(x-.65).^2/6) - (x>.225&x<.275).*exp(-(x-.25).^2/6));
end

function chi = horizon_function(x,ell,xi)
  chi = zeros(size(xi));
  
  for i=1:length(chi)
    chi = ( xi(i)>=x-ell ) && ( xi(i)<=x+ell );
  end
end
    
% function f = f_function(eta,xi)
%   f = sign(eta+xi).*f_tilde(abs(eta+xi),xi);
% end
% function ft = f_tilde(s,xi)
%   c     = 0.3;
%   delta = 0.2;
%   ft    = c*s.*(norm(xi)<delta);
% end

function f = f_function(eta,xi)
  y = abs(eta+xi);
  s = (y-abs(xi))./abs(xi+0.00001);
  k = 1;
  delta = .3;
  c = 18*k/(pi*delta^4);
  f = c*s;

  f = sign(eta+xi).*f;
%   delta = .3;
%   E     = 1e3;
%   c     = 4*E/delta^3/sqrt(pi);
%   f     = c*exp(-xi.^2/delta^2).*eta;
end


%  Also requires the following functions from 'fem_functions':
%
%  oned_bilinear
%  oned_f_int
%  oned_gauss
%  oned_mesh
%  oned_shape
