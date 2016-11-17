function [nu, freq_cm, num_onda, num_onda_cm, perfil, s, nome_fig, M, M_from, M_to, lambda, elem, lim_perfil] = main(h161, h163)
   %%----------------------------------------------------------------
   %% Dados de entrada:
   elem = "Dy";
#   lambda = 599.022; %% lambda vacuo que é 598.9 no lambda ar / 16693.87 cm-1
   lambda = 597.614; %% lambda vacuo que é 597.452 no lambda ar / 16733 cm-1

   lim_perfil = 0.001;
   nome_fig = ["plot_" elem "_simulado_" num2str(lambda)];

   v_luz = 299792458; %% velocidade da luz (m/s) no vacuo


   c_cons = 1e-8; %1.5; %0.00000001; %1.5; % peso do pedestal gaussiano
   gama = 20; %1; %20; %% largura da lorentziana à meia altura em MHz (FWHM = 2 * HWHM)

   T = 1300; %% 1300.0 / (4*log(2)); %% temperatura em Kelvin (influencia no pedestal gaussiano)
   M = [156.0, 158.0, 160.0, 162.0, 164.0]; %% massas dos isotopos

   M_from = [0.0,0.0,0.0,0.0,0.0];
   M_to = [0.0,0.0,0.0,0.0,0.0];

   %% s0 são as abundancias naturais dos isótopos
   s0 = [
   0.00056, %% 156
   0.00095, %% 158
   0.02329, %% 160
   0.18889, %% 162
   0.25475 %% 164
   ];

   s0_impar = [
   0.24896, %% 161
   0.2826 %% 163
   ];



#   h161(:,2)
#   h163(:,2)
#   sum(h161(:,2))
#   sum(h163(:,2))
   


   #delta_nu = 1000;
   
   h161(:,2) = h161(:,2) .* s0_impar(1)/sum(h161(:,2));
   h163(:,2) = h163(:,2) .* s0_impar(2)/sum(h163(:,2));
   
   
   h161(:,2)
   h163(:,2)
   sum(h161(:,2))
   sum(h163(:,2))
   
   
   
   
   # ajustando os vetores s0 e M para comportar a parte da estrutura hiperfina.
   s0 = [s0; h161(:,2); h163(:,2)];
   M = [M, 161*ones(1, length(h161(:,2))), 163*ones(1, length(h163(:,2)))];
   
   M_from = [M_from, h161(:,3)', h163(:,3)'];
   M_to = [M_to, h161(:,4)', h163(:,4)'];
   
#   %% shifts isotopicos (MHz) para lambda = 599,022 nm no vacuo = 598,9 nm no ar (16693 cm-1)
#   sh = [
#   -1359.6, %% 156-158
#   -923.6, %% 158-160
#   -932.8, %% 160-162
#   -866.9 %% 162-164
#   ];
#   
#   
##   delta_impar_163 = -1000;
##   delta_impar_161 = -500;
#   delta_impar_163 = 0;
#   delta_impar_161 = 0;
#   %% shifts isotopicos (MHz) para lambda = 599,022 nm no vacuo = 598,9 nm no ar (16693 cm-1)
#   sh_impar = [
#   -251.4 + delta_impar_161, %% 160-161
#   -282.6 + delta_impar_163 %% 162-163
#   ];


   %% shifts isotopicos (MHz) para lambda = 597,614 nm no vacuo = 597,4 nm no ar (16733 cm-1)
   sh = [
   -1547.1, %% 156-158
   -1037.9, %% 158-160
   -1055.0, %% 160-162
   -977.2 %% 162-164
   ];

   delta_impar_163 = 0;
   delta_impar_161 = 0;
   %% shifts isotopicos (MHz) para lambda = 597,614 nm no vacuo = 597,4 nm no ar (16733 cm-1)
   %% Ref: Wakasugi 1990
   sh_impar = [
   -278.6 + delta_impar_161, %% 160-161
   -311.5 + delta_impar_163 %% 162-163
   ];
   %% shifts isotopicos (MHz) para lambda = 597,614 nm no vacuo = 597,4 nm no ar (16733 cm-1)
   %% Ref: tese do Eliel 1979
#   sh_impar = [
#   780.0 + delta_impar_161, %% 161-162
#   -308.0 + delta_impar_163 %% 162-163
#   ];






   %% shifts isotopicos acumulados
   sh_ac(1) = sh(1); %% 158
   sh_ac(2) = sum(sh(1:2)); %% 160
   sh_ac(3) = sum(sh(1:3)); %% 162
   sh_ac(4) = sum(sh(1:4)); %% 164
   
   sh_impar(1) = sh_ac(2) + sh_impar(1); %% 161 - Se for usar Wakasugi
#   sh_impar(1) = sh_ac(3) + sh_impar(1); %% 161 - Se for usar Eliel / Zaal
   sh_impar(2) = sh_ac(3) + sh_impar(2); %% 163

# Essa concatenação tem que estar abaixo dos cálculos dos shifts dos impares.
# O erro que tivemos durante mais de 1 ano foi devido a isso.
   sh_ac = [0, sh_ac];
   
      
   ## fórmula para cálculo do centro de massa: cm = sum(perfil .* nu) / sum(perfil)
   ## usamos essa fórmula abaixo:
#   cm161 = sum(h161(:,2) .* (h161(:,1).-h161(1,1))) / sum(h161(:,2)); %% centro de massa da transição h161
#   cm163 = sum(h163(:,2) .* (h163(:,1).-h163(1,1))) / sum(h163(:,2)); %% centro de massa da transição h163
   
   
##  Mudanças relativas ao cálculo do centro de massa:

#   h161_norm = h161(:,2) / max(h161(:,2));
#   h163_norm = h163(:,2) / max(h163(:,2));
#   
#   # Cálculos do centro de massa somente para as linhas mais relevantes.
#   cm161 = sum(h161(h161_norm >= 0.2, 2) .* h161(h161_norm >= 0.2, 1)) / sum(h161(h161_norm >= 0.2, 2)); %% centro de massa da transição h161
#   cm163 = sum(h163(h163_norm >= 0.2, 2) .* h163(h163_norm >= 0.2, 1)) / sum(h163(h163_norm >= 0.2, 2)); %% centro de massa da transição h163

#   disp("size h161");
#   disp(size(h161(h161_norm >= 0.1, 2)));
#   
#   disp("size h163");
#   disp(size(h163(h163_norm >= 0.1, 2)));


# Cálculo completo dos centros de massa:   
#   cm161 = sum(h161(:,2) .* h161(:,1)) / sum(h161(:,2)); %% centro de massa da transição h161
#   cm163 = sum(h163(:,2) .* h163(:,1)) / sum(h163(:,2)); %% centro de massa da transição h163
   
   cm161 = 0;
   cm163 = 0;
   
#   disp("size h161");
#   disp(size(h161));
#   
#   disp("size h163");
#   disp(size(h163));
   
   
#   disp("cm161:");
#   printf("%.15f\n", cm161);
#   disp("cm163");
#   printf("%.15f\n", cm163);
   
   h161(:,1) = h161(:,1) + sh_impar(1) - cm161;
   h163(:,1) = h163(:,1) + sh_impar(2) - cm163; 

#   h161(:,1) = h161(:,1) + sh_impar(1);
#   h163(:,1) = h163(:,1) + sh_impar(2);
   
#   disp("sh_ac:");
#   disp(sh_ac);
   
   sh_ac = [sh_ac'; h161(:,1); h163(:,1)];
   
   
   [nu, freq_cm, num_onda, num_onda_cm, perfil, s, nome_fig, M, lambda, elem, lim_perfil, cm] = shifts(elem, lambda, lim_perfil, c_cons, gama, T, M, sh_ac, s0);
   
#   cm161
#   cm163
   
#   disp("Shift 161:");
#   disp(sh_impar(1));
#   disp("Shift 163:");
#   disp(sh_impar(2));
#   
#   disp(["freq_cm = " num2str(freq_cm)]);
#   disp(["CM = " num2str(cm)]);
   
   cm161 = sh_impar(1) + freq_cm - cm;
   cm163 = sh_impar(2) + freq_cm - cm;
   
#   disp("cm161:");
#   printf("%.15f\n", cm161);
#   disp("cm163");
#   printf("%.15f\n", cm163);
   
   
endfunction
