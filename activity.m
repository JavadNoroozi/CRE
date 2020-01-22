
function gamma_i = activity(m_i)
  
    global A
   global Z 
    
    I = 0.5*sum(m_i.*Z.^2);
    gamma_i = exp(-A*Z.^2*log(10)*(sqrt(I)/(1 + sqrt(I))-0.3*I));
    c = (A*log(10)/55.508)*(2*(I + 2*sqrt(I))/(1 + sqrt(I)) - 4*log(1 + sqrt(I))-0.3*I^2);
    a_i = m_i.*gamma_i;

    #  water activity mu(H2O)=mu*+ RTln aw = mu(dagger)+ RTln(55.34)+RT ln aw
    #  = mu(dagger) + RT ln (aw*55.34)
    # We may call aw*55.34 apparent activity of the water given by
    # ln aw = ln (55.34) + c-(1-xw)/xw
    a_i(1) = exp( c + 1 - sum(m_i)/m_i(1));
    gamma_i = a_i./m_i;
 end
