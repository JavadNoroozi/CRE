function constrains = mass_charge(m_vector)
  global K1
  global K2
  global K3
  global K4
  global K5
  global M
  global alpha
  
  
# H2O:1 RNH2:2 CO2:3 RNH3:4 RNHCOO:5 HCO3:6 H3O:7 OH:8 CO3(-2):9

  m_H2O = m_vector(1);
  m_RNH2 = m_vector(2);
  m_CO2 = m_vector(3);
  m_RNH3 = m_vector(4);
  m_RNHCOO = m_vector(5);
  m_HCO3 = m_vector(6);
  m_H3O = m_vector(7);
  m_OH = m_vector(8);
  m_CO3 = m_vector(9);
  
  
  gamma_vector = activity(m_vector);
  #gamma_vector = [ 1, 1, 1, 1, 1, 1, 1, 1, 1] ;
  
  gamma_H2O = gamma_vector(1);
  gamma_RNH2 = gamma_vector(2);
  gamma_CO2 = gamma_vector(3);
  gamma_RNH3 = gamma_vector(4);
  gamma_RNHCOO = gamma_vector(5);
  gamma_HCO3 = gamma_vector(6);
  gamma_H3O = gamma_vector(7);
  gamma_OH = gamma_vector(8);
  gamma_CO3 = gamma_vector(9);
  

constrains(1) = M -(m_RNH2 + m_RNHCOO + m_RNH3); # amine balance
constrains(2) = M*alpha -( m_RNHCOO + m_HCO3 + m_CO2 + m_CO3); # CO2 balance
constrains(3) = (m_H3O + m_RNH3) - (m_RNHCOO + m_HCO3 + 2*m_CO3 + m_OH); # chage balance

# Reaction equilibrium
# R1: RNH3 + RNHCOO = 2 RNH2 + CO2
# R2: RNH3 + H2O = RNH2 + H3O+
# R3: RNHCOO + H2O = RNH2 + HCO3
# R4: RNH2 + H2O = RNH3 + OH
# R5 : HCO3(-) + H2O = CO3(-2) + H3O(+)


constrains(4) = exp(-K1)- (m_RNH2*m_RNH2*m_CO2)/(m_RNH3*m_RNHCOO)  *  (gamma_RNH2*gamma_RNH2*gamma_CO2)/(gamma_RNH3*gamma_RNHCOO) ;
constrains(5) = exp(-K2)- (m_RNH2*m_H3O)/(m_RNH3) *  (gamma_RNH2*gamma_H3O)/(gamma_RNH3);
constrains(6) = exp(-K3) - (m_RNH2*m_HCO3)/(m_RNHCOO) * (gamma_RNH2*gamma_HCO3)/(gamma_RNHCOO);
constrains(7) = exp(-K4) - (m_OH*m_RNH3)/(m_RNH2) * (gamma_OH*m_RNH3)/(gamma_RNH2);
constrains(8) = exp(-K5) - (m_H3O*m_CO3)/(m_HCO3) * (gamma_H3O*gamma_CO3)/(gamma_HCO3);


end
