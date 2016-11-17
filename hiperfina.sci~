// Esse programa fornece os valores das possíveis transições da estrutura hiperfina, bem como calcula suas intensidades relativas.
// A explicação dos conceitos utilizados nos cálculos encontram-se em:

// Kuhn, H. G. Atomic Spectra. Cap. 6 - Hyperfine Structure, pgs. 191 e 333, Longmans, 1969. 

// Destro, M. G. Espectroscopia a Laser em Vapor Metálico de Urânio. Cap. 2, pg. 13, Tese - ITA.

// Harrison - Practical Spectroscopy - pg. 248

//-------------------------------------------------------------
function W = calcW(A, B, cond, F, I, J)

   // Equações de Casimir:
	// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
	C = F(cond).*(F(cond)+1) - I.*(I+1) - J.*(J+1);
   
   // Equação de Casimir
   W = A .* C./2 + B .* (3.*C.*(C+1) - 4 .*I.*(I+1).*J.*(J+1)) ./ (8 .* I .*(2.*I-1).* J.*(2.*J-1));
		
endfunction
//-------------------------------------------------------------

function [resultado] = hiperfina(I, J, J_sup, Ai, Af, Bi, Bf)
//----------------------------------------------
// Coloque aqui os valores de I, J (J de partida da transição) e J_sup, que é o J de destino da transição.

//	I = 7/2;
//	J = 6;
//	J_sup = 7;

	//----------------------------------------------
	v1 = 0:2*I;

//	F = abs(v1 + J - I);
//	F_sup = abs(v1 + J_sup - I);

//	F = (v1 + J - I);
//	F_sup = (v1 + J_sup - I);

	F = abs(J - I) : abs(J+I);
	F_sup = abs(J_sup - I) : abs(J_sup+I);

	[nl1, nc1] = size(F);
	[nl2, nc2] = size(F_sup);
	//disp(size(F));


	matF = F';
	for i=1:nc2-1
		matF = [matF, F'];
	end
//	disp(matF);
//	disp(size(matF));
	matF_sup = F_sup;
	for i=1:nc1-1
		matF_sup = [matF_sup; F_sup];
	end
//	disp(matF_sup);
//	disp(size(matF_sup));
//	disp(size(F));
//	disp(F);
//	disp(size(F_sup));
//	disp(F_sup);
//	abort;
	transicoes = matF_sup - matF;
	//disp(transicoes);

	criterio = transicoes == 0 | transicoes == 1 | transicoes == -1;
	//disp(criterio);

	partida = matF(criterio);
	chegada = matF_sup(criterio);

	//disp(partida);
	//disp(chegada);

	//disp("partida, chegada:");
	//disp([partida, chegada]');


	if J_sup-J == 1
		// IF (intensidade da transicao J-1 -> J) obs: usa-se J_sup nesse caso (o J do nível superior).
		// Fórmula da pág. 13 da tese do Marcelo.
		
		//-----------------------------------------------------------
		// transicao F-1 -> F
		cond = partida == chegada-1;
		//disp("cond:");
		//disp(cond');
		IF_menos1_F = (J_sup+chegada(cond)+I+1) .* (J_sup+chegada(cond)+I) .* (J_sup+chegada(cond)-I) .* (J_sup+chegada(cond)-I-1) ./ chegada(cond);

      // Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW_menos1 = Wf-Wi;
      //-----------------------------------------------------------
      
		// transicao F -> F
		cond = partida == chegada;
		//disp("cond:");
		//disp(cond');
		IF_F = -(J_sup+chegada(cond)+I+1) .* (J_sup+chegada(cond)-I) .*(J_sup-chegada(cond)+I) .* (J_sup-chegada(cond)-I-1) .* (2*chegada(cond)+1) ./(chegada(cond).^2 + chegada(cond));

      // Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW = Wf-Wi;
      //-----------------------------------------------------------
      
		// transicao F+1 -> F
		cond = partida == chegada+1;
		//disp("cond:");
		//disp(cond');
		IF_mais1_F = (J_sup-chegada(cond)+I) .* (J_sup-chegada(cond)+I-1) .* (J_sup-chegada(cond)-I-1) .* (J_sup-chegada(cond)-I-2) ./ (chegada(cond)+1);
		//disp("J -> J+1");
		
		// Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW_mais1 = Wf-Wi;
		//-----------------------------------------------------------
	elseif J_sup == J

		// IF (intensidade da transicao J -> J) obs: usa-se J_sup nesse caso (o J do nível superior).
		// Fórmula da pág. 13 da tese do Marcelo.
		
		//-----------------------------------------------------------
		// transicao F-1 -> F
		cond = partida == chegada-1;
		//disp("cond:");
		//disp(cond');
		IF_menos1_F = -(J_sup+chegada(cond)+I+1) .* (J_sup+chegada(cond)-I) .* (J_sup-chegada(cond)+I+1) .* (J_sup-chegada(cond)-I) ./ chegada(cond);
		//IF_menos1_F = -(J_sup+partida(cond)+I+1) .* (J_sup+partida(cond)-I) .* (J_sup-partida(cond)+I+1) .* (J_sup-parida(cond)-I) ./ partida(cond);
      
      // Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW_menos1 = Wf-Wi;
      //-----------------------------------------------------------
      
		// transicao F -> F
		cond = partida == chegada;
		//disp("cond:");
		//disp(cond');
		IF_F = (J_sup .* (J_sup+1) + chegada(cond) .* (chegada(cond)+1) - I .* (I+1)).^2 .* (2 * chegada(cond) + 1) ./ (chegada(cond) .* (chegada(cond)+1));

      // Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW = Wf-Wi;
      
		//-----------------------------------------------------------
		// transicao F+1 -> F
		cond = partida == chegada+1;
		//disp("cond:");
		//disp(cond');
		IF_mais1_F = -(J_sup+chegada(cond)+I+2) .* (J_sup+chegada(cond)-I+1) .* (J_sup-chegada(cond)+I) .* (J_sup-chegada(cond)-I-1) ./ (chegada(cond)+1);
		//IF_mais1_F = -(J_sup+partida(cond)+I+2) .* (J_sup+partida(cond)-I+1) .* (J_sup-partida(cond)+I) .* (J_sup-partida(cond)-I-1) ./ (partida(cond)+1);
		//disp("J -> J");
		
		// Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW_mais1 = Wf-Wi;
		//-----------------------------------------------------------
		
	elseif J_sup - J == -1
	
		// IF (intensidade da transicao J -> J-1) obs: usa-se J_sup nesse caso como o J do nível superior. Veja que nesse caso as fórmulas usam J ao invés de J_sup.
		// Fórmula da pág. 13 da tese do Marcelo e do livro do Harrison-Practical Spectroscopy pag.248/268 no djvu.
		
		//-----------------------------------------------------------
		// transicao F-1 -> F
		cond = partida == chegada+1;
		//disp("cond:");
		//disp(cond');
		IF_mais1_F = (J+partida(cond)+I+1) .* (J+partida(cond)+I) .* (J+partida(cond)-I) .* (J+partida(cond)-I-1) ./ partida(cond);
      
      // Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW_mais1 = Wf-Wi;
      
      //-----------------------------------------------------------
		// transicao F -> F
		cond = partida == chegada;
		//disp("cond:");
		//disp(cond');
		IF_F = -(J+partida(cond)+I+1) .* (J+partida(cond)-I) .*(J-partida(cond)+I) .* (J-partida(cond)-I-1) .* (2*partida(cond)+1) ./(partida(cond).^2 + partida(cond));
		
		// Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW = Wf-Wi;
      
      //-----------------------------------------------------------
		// transicao F+1 -> F
		cond = partida == chegada-1;
		//disp("cond:");
		//disp(cond');
		IF_menos1_F = (J-partida(cond)+I) .* (J-partida(cond)+I-1) .* (J-partida(cond)-I-1) .* (J-partida(cond)-I-2) ./ (partida(cond)+1);
		
		// Equações de Casimir:
		// Equação da teoria (Casimir) é: C = F(f+1) - I(I+1) - J(J+1)
		clear Wi;
		clear Wf;
      Wi = calcW(Ai, Bi, cond, partida, I, J);
      Wf = calcW(Af, Bf, cond, chegada, I, J_sup);
      
      deltaW_menos1 = Wf-Wi;
		//-----------------------------------------------------------
	end

	maximo = max([IF_menos1_F; IF_F; IF_mais1_F]);


//	disp([IF_menos1_F*100/maximo, partida(partida == chegada-1), chegada(partida == chegada-1)]);
//	disp([IF_F*100/maximo, partida(partida == chegada), chegada(partida == chegada)]);
//	disp([IF_mais1_F*100/maximo, partida(partida == chegada+1), chegada(partida == chegada+1)]);
	
	IF_menos1_F_norm = IF_menos1_F/maximo;
	IF_F_norm = IF_F/maximo;
	IF_mais1_F_norm = IF_mais1_F/maximo;
	
	
	ret1 = [IF_menos1_F_norm, partida(partida == chegada-1), chegada(partida == chegada-1)];
	
	ret2 = [IF_F_norm, partida(partida == chegada), chegada(partida == chegada)];
	
	ret3 = [IF_mais1_F_norm, partida(partida == chegada+1), chegada(partida == chegada+1)];
	
	resultado = [[deltaW_menos1; deltaW; deltaW_mais1], [ret1; ret2; ret3]];
	
//	disp([resultado]);
//	disp(deltaW_menos1);
//	disp(deltaW);
//	disp(deltaW_mais1);
//	disp([[IF_menos1_F_norm; IF_F_norm; IF_mais1_F_norm], partida, chegada]);
//------------------------------------------------------------	
	// montagem da matriz F versus F_sup com os valores dessas transicoes
//	disp(resultado);
	
//	[nl, nc] = size(criterio);
////	disp([nl, nc]);
//	//criterio = transicoes == 0 | transicoes == 1 | transicoes == -1;

////	disp(criterio);
//	

//	cont1 = 1;
//	cont2 = 1;
//	cont3 = 1;
//	
//	for i=1:nl
//		for j=1:nc
//		
//			if(F_sup(j) == F(i))
//				matriz_final(i, j) = ret2(cont2, 1);
//				cont2 = cont2+1;
//			
//			elseif(F_sup(j)-1 == F(i))
//				matriz_final(i, j) = ret1(cont1);
//				cont1 = cont1+1;
//				
//			elseif(F_sup(j)+1 == F(i))
//				matriz_final(i, j) = ret3(cont3);
//				cont3 = cont3+1;
//			end
//			
//		end
//	end 
//	
//	matriz_final = [F_sup; matriz_final];
//	matriz_final = [[0; F'], matriz_final];
//	
//	matriz_final = string(matriz_final);
//	matriz_final(matriz_final == '0') = '';
//	vetJ_sup(1:nl+1) = '';
////	disp(vetJ_sup);
//	vetJ_sup(1) = "J = " + string(J_sup)+ '; F''';
//	[nl, nc] = size(matriz_final);
//	matriz_final = [matriz_final(:, 1), vetJ_sup, matriz_final(:, 2:nc)];
//	
//	[nl, nc] = size(matriz_final);
////	disp(nl, nc);
////	disp(size(matriz_final));
//	vet_vazio(1:nc) = '';
////	disp(vet_vazio');
//	matriz_final = [matriz_final(1, :); vet_vazio'; matriz_final(2:nl, :)];
//	
//	matriz_final(2, 1) = "J = " + string(J)+ '; F';
//	//disp(matriz_final);
//	
//	[nl, nc] = size(matriz_final);
////	aqui redimensiona-se a matriz para não dar erro nos próximos J. 
//	clear vet_vazio;
//	if nl > nc
//		vet_vazio(1:nl, 1:nl-nc) = '';
//		matriz_final = [matriz_final, vet_vazio];
//	elseif nc > nl
//		vet_vazio(1:nc-nl, 1:nc) = '';
//		matriz_final = [matriz_final; vet_vazio];
//	end
	
//----------------------------------------------------------- 

endfunction
