function [K,f] = elasticity_matrices(x,e_conn,...
                                     material     )
%-------------------------------------------------------------------------------
%  Computes the two-dimensional finite element matrices for the Lame
%  problem.
  young = material.young;
  nu    = material.nu;
  
  lam   = nu*young/((1+nu)*(1-2*nu)); 
  mu    = young/(1+nu)/2;
  
  n_nodes              = size(x,1); 
  [n_elements,nel_dof] = size(e_conn);

  [rr,ss,wt]  = twod_gauss(7);
  one         = ones(size(wt));
  n_equations = 2*n_nodes;                  % build everything and deflate later
  ide_u = [(1:2:2*n_nodes-1)' 2*(1:n_nodes)'];
  
  % preallocating space
  II        = zeros((2*nel_dof)^2*n_elements,1);
  JJ        = zeros((2*nel_dof)^2*n_elements,1);
  XX        = zeros((2*nel_dof)^2*n_elements,1);
  
  f         = zeros (n_equations, 1);       % weak residual vector contributions
  
  entry_counter = 0;
  for n_el=1:n_elements
    % compute value of each test function and their spatial derivaties
    % at the integration points.
    nodes_local                = e_conn(n_el,:);
    x_local                    = x(nodes_local,:);
    [x_g,wt_g,phi,phi_x,phi_y] = twod_shape(x_local,rr,ss,wt);

    [fx,fy] = force(x_g);
    
    k11_loc = (2*mu+lam)*twod_bilinear( one, phi_x, phi_x, wt_g) ...
            +    mu     *twod_bilinear( one, phi_y, phi_y, wt_g);  
    k12_loc =    mu     *twod_bilinear( one, phi_x, phi_y, wt_g) ...
            +       lam *twod_bilinear( one, phi_y, phi_x, wt_g);
    k21_loc =    mu     *twod_bilinear( one, phi_y, phi_x, wt_g) ...
            +       lam *twod_bilinear( one, phi_x, phi_y, wt_g);
    k22_loc =    mu     *twod_bilinear( one, phi_x, phi_x, wt_g) ...
            + (2*mu+lam)*twod_bilinear( one, phi_y, phi_y, wt_g);
    
    f1_loc  = twod_f_int( fx, phi, wt_g );
    f2_loc  = twod_f_int( fy, phi, wt_g );
    
    for n_t=1:nel_dof
      n_test = ide_u(nodes_local(n_t),1);
      for n_u=1:nel_dof
        n_unku = ide_u(nodes_local(n_u),1);
        n_unkv = ide_u(nodes_local(n_u),2);

        entry_counter = entry_counter + 1;
        II( entry_counter ) = n_test;
        JJ( entry_counter ) = n_unku;
        XX( entry_counter ) = k11_loc(n_t,n_u);
        %
        entry_counter = entry_counter + 1;
        II( entry_counter ) = n_test;
        JJ( entry_counter ) = n_unkv;
        XX( entry_counter ) = k12_loc(n_t,n_u);
      end
      f(n_test) = f(n_test) + f1_loc(n_t);
      
      n_test = ide_u(nodes_local(n_t),2);
      for n_u=1:nel_dof
        n_unku = ide_u(nodes_local(n_u),1);
        n_unkv = ide_u(nodes_local(n_u),2);

        entry_counter = entry_counter + 1;
        II( entry_counter ) = n_test;
        JJ( entry_counter ) = n_unku;
        XX( entry_counter ) = k21_loc(n_t,n_u);
        %
        entry_counter = entry_counter + 1;
        II( entry_counter ) = n_test;
        JJ( entry_counter ) = n_unkv;
        XX( entry_counter ) = k22_loc(n_t,n_u);
      end
      
      f(n_test) = f(n_test) + f2_loc(n_t);
    end
  end % element loop
 
  K = sparse( II(1:entry_counter), JJ(1:entry_counter), ...
              XX(1:entry_counter),                      ...
              n_equations,         n_equations             );

end % function element_integrations

function [fx,fy] = force(x)
  n  = size(x,1);
  fx = zeros(n,1);
  fy = zeros(n,1);
end