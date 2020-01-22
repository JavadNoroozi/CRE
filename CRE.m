# initial guess
# m_CO2,m_RNH2, m_RNH3,m_RNHCOO, m_HCO3, m_H3O,m_OH,m_CO3
clc
clear all
global alpha
global  M
global K1 
global K2 
global K3 
global K4
global K5
global A
global Z


A = 0.35;  
Z = [0, 0, 0, 1, -1, -1, 1, -1, -2] ; # species charge valence

    
WF = 0.3 ; # weight fraction o amine
MW_H2O = 18.01528; # molecular weight of water and amine
MW_amine = 61.08; # molecular weight of amine
M = ((WF/(1-WF))/MW_amine)*1000;


Delta_pka = 0;
Delta_pka_carb = 0;
 # R1: RNH3 + RNHCOO = 2 RNH2 + CO2
# R2: RNH3 + H2O = RNH2 + H3O+
# R3: RNHCOO + H2O = RNH2 + HCO3
# R4: RNH2 + H2O = RNH3 + OH
# R5 : HCO3(-) + H2O = CO3(-2) + H3O(+)

K1 = 14.87 + Delta_pka + Delta_pka_carb ;         # combining Pka, carbamate formation and HCO3-
K2 = 21.88 + Delta_pka ;        # PKa reaction
K3 = 3.65 + Delta_pka_carb  ;  # Carbmate reaction
K4 = 10.38 - Delta_pka   ;  # Water ionization and Pka
K5 = 23.80  ;    # CO3(-2)



m_initial_guess = [5.5508e+01   6.3144e+00   6.9939e-10   3.6120e-01   3.4091e-01   1.4033e-03   1.7991e-11   1.8638e-03   8.5108e-03];

alpha = 0.05;

for i=1:100

alpha = alpha + 0.01;

#m_initial_guess = [5.55084351e+01,6.25268924e+00,8.12950410e-10,3.83407992e-01, 3.66180998e-01, 1.34212056e-03, 1.70046803e-11, 6.87015708e-04,7.59892852e-03];

options = optimset ('TolX', 10E-30,'TolFun', 10E-30,'MaxIter',1000000000000 );
 [x_found,zero] = fsolve(@mass_charge , m_initial_guess, options);
 molality(i,:) = x_found;
 mole_fractions(i,:)= molality(i,:)/sum(molality(i,:));
 m_initial_guess = x_found;
 loading(i)=alpha;

end


Co2_pressure = mole_fractions(:,3).* 200*10E3;

#plot(loading,mole_fractions(:,2))
#hold on
#plot(loading,mole_fractions(:,3))
#hold on
#plot(loading,mole_fractions(:,4))
#hold on
#plot(loading,mole_fractions(:,5))
#hold on
#plot(loading,mole_fractions(:,6))
#hold on
#plot(loading,mole_fractions(:,7))
#hold on
#plot(loading,mole_fractions(:,8))
#hold on
#plot(loading,mole_fractions(:,9))
#hold on
semilogy(loading,Co2_pressure)
hold on








