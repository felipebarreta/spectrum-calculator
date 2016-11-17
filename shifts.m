#Relações de gaussianas e lorentzianas (e suas equações):

#http://mathworld.wolfram.com/GaussianFunction.html

#http://mathworld.wolfram.com/LorentzianFunction.html

#http://mathworld.wolfram.com/FullWidthatHalfMaximum.html


function [nu, freq_cm, num_onda, num_onda_cm, perfil, s, nome_fig, M, lambda, elem, lim_perfil, cm] = shifts(elem, lambda, lim_perfil, c_cons, gama, T, M, sh_ac, s0)
   %%----------------------------------------------------------------
   %% Dados de entrada:
#   elem = "Dy";
#   lambda = 599.022; %% lambda vacuo que é 598.9 no lambda ar
#   #lambda = 597.614; %% lambda vacuo que é 597.94 no lambda ar

#   lim_perfil = 0.001;
   nome_fig = ["plot_" elem "_simulado_" num2str(lambda)];

   v_luz = 299792458; %% velocidade da luz (m/s) no vacuo


#   c_cons = 1.5;
#   gama = 70; %% largura da lorentziana à meia altura em MHz (FWHM = 2 * HWHM)

#   T = 1500; %% temperatura em Kelvin
#   M = [156, 158, 160, 162, 164, 161, 163]; %% massas dos isotopo


   #delta_nu = 1000;


   %% shifts isotopicos (MHz) para lambda = 599,022 nm no vacuo = 598,9 nm no ar
#   sh = [
#   -1359.616, %% 156-158
#   -923.613, %% 158-160
#   -932.815, %% 160-162
#   -866.912 %% 162-164
#   ];
#   %% shifts isotopicos (MHz) para lambda = 599,022 nm no vacuo = 598,9 nm no ar
#   sh_impar = [
#   -251.418, %% 160-161
#   -282.62 %% 162-163
#   ];

   #%% shifts isotopicos (MHz) para lambda = 597,614 nm no vacuo = 597,4 nm no ar
   #sh = [
   #-1547.124, %% 156-158
   #-1037.924, %% 158-160
   #-1055.02, %% 160-162
   #-977.214 %% 162-164
   #];

   #%% shifts isotopicos (MHz) para lambda = 597,614 nm no vacuo = 597,4 nm no ar
   #sh_impar = [
   #-278.622, %% 160-161
   #-311.521 %% 162-163
   #];

   %% s0 são as abundancias naturais dos isótopos
#   s0 = [
#   0.00056, %% 156
#   0.00095, %% 158
#   0.02329, %% 160
#   0.18889, %% 162
#   0.25475 %% 164
#   0.24896/10 %% 161
#   0.2826/10 %% 163
#   ];

   %%----------------------------------------------------------------
   freq_cm = 1000 * v_luz / lambda; %% frequência da onda em MHz do centro de massa da transição
   delta_nu = 2 * 3.581e-7 * freq_cm * sqrt(T./M);

   minimo = min(sh_ac)-1000*gama;
   maximo = max(sh_ac)+1000*gama;
   nu = linspace(minimo, maximo, 10000);


   length(sh_ac)
   s = [];
   gama_hwhm = gama/2;
   for i = 1:length(sh_ac)
	   nu0 = sh_ac(i);

# Equação para Alargamento Doppler FWHM 
	   s_temp = s0(i) .* ((gama^2 / 4) ./ ((nu - nu0).^2 + gama^2 / 4) + c_cons .* exp(-4*log(2) * (nu - nu0).^2 ./ delta_nu(i)^2)) .* exp(-4*log(2) * (nu - nu0).^2 ./ delta_nu(i)^2);

# Equação para Alargamento Doppler full width at 1/e maximum
#s_temp = s0(i) .* ((gama^2 / 4) ./ ((nu - nu0).^2 + gama^2 / 4) + c_cons .* exp(-(nu - nu0).^2 ./ delta_nu(i)^2)) .* exp(-(nu - nu0).^2 ./ delta_nu(i)^2);

# Equação de acordo com Hansch e Smith, mas com sinal neativo na exponencial
#      s_temp = s0(i) .* ((gama_hwhm^2) ./ (4*pi^2*(nu - nu0).^2 + gama_hwhm^2) + c_cons .* exp(-(nu - nu0).^2 ./ delta_nu(i)^2)) .* exp(-(nu - nu0).^2 ./ delta_nu(i)^2);
	   s = [s; s_temp];	
   endfor


   [nl,nc] = size(s);

   perfil = zeros(1, nc);
   for i = 1:nl
	   perfil = perfil + s(i,:);
   endfor


   s = s./ max(perfil); %% normalização de s

   perfil = perfil ./ max(perfil); %% normalização do perfil

   cm = sum(perfil .* nu) / sum(perfil); %% centro de massa da transição

   nu = nu + freq_cm - cm; %% calculo considerando o desvio do centro de massa e o centro da transição
   
   num_onda = nu * 1e4 / v_luz;
   num_onda_cm = freq_cm * 1e4 / v_luz;

#   callpy("plot", "plot", nu, freq_cm, num_onda, num_onda_cm, perfil, s, nome_fig, M, lambda, elem, lim_perfil);

   #hold on;

   #plot(num_onda, s(1,:));
   #plot(num_onda, s(2,:));
   #plot(num_onda, s(3,:));
   #plot(num_onda, s(4,:));
   #plot(num_onda, perfil);
   #plot(num_onda_cm, 0.5);

   #hold off;
endfunction
