// A única coisa que essa função faz é chamar a função hiperfina.sci. Foi um jeito de
// contornar o problema ValueError: too many values to unpack do scilab2py, uma vez
// que no arquivo hiperfina.sci há duas funções definidas.

function [resultado] = callHiperfina(I, J, J_sup, Ai, Af, Bi, Bf)

   // Esse funcprot só evita um warning acerca da redefinição da função hiperfina.sci.
   funcprot(0);

   // Obs: Esse arquivo hiperfina.sci é o arquivo que resolve todos os casos que queremos, ou seja, deltaJ = +1, 0, -1.
   exec('hiperfina.sci');
   
   disp resultado;
   // chamada da função hiperfina.sci
   resultado = hiperfina(I, J, J_sup, Ai, Af, Bi, Bf);
   
endfunction
