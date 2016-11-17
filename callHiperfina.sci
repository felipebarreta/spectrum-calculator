//# <one line to give the program's name and a brief idea of what it does.>
//#    Copyright (C) 2016  Luiz Felipe Nardin Barreta

//#    This program is free software: you can redistribute it and/or modify
//#    it under the terms of the GNU General Public License as published by
//#    the Free Software Foundation, either version 3 of the License, or
//#    (at your option) any later version.

//#    This program is distributed in the hope that it will be useful,
//#    but WITHOUT ANY WARRANTY; without even the implied warranty of
//#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//#    GNU General Public License for more details.

//#    You should have received a copy of the GNU General Public License
//#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



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
