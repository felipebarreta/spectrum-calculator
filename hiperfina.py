#!/usr/bin/python
# -*- coding: utf-8 -*-


#---------------------------------------------------------------------
# imports:
from scilab2py import scilab as sci
from oct2py import octave as oc
import numpy as np
import matplotlib.pyplot as plt
import os
import logging

#---------------------------------------------------------------------
# Variaveis globais (inicialização):
I = 0.0
J = 0.0
J_sup = 0.0
Ai = 0.0
Bi = 0.0
Af = 0.0
Bf = 0.0
filename = ""
#---------------------------------------------------------------------
def plot(nu, freq_cm, num_onda, num_onda_cm, perfil, s, nome_fig, M, M_from, M_to, lamb, elem, lim_perfil):

   
   nl, nc = s.shape
   
   picosx = []
   picosy=[]
   picosy2=[]
   for i in range(nl):
      picosx.append(num_onda[0][np.where(s[i,:] == s[i,:].max())])
      picosy.append(perfil[0][s[i,:] == s[i,:].max()])
      picosy2.append(s[i][np.where(s[i,:] == s[i,:].max())])
      
#   print "picosx = "
#   print picosx
#   
#   print "picosy = "
#   print picosy

   plt.rcParams['font.family'] = 'Times New Roman' # aqui a gente seta que vai usar fonte times new roman no gráfico
   fig = plt.figure()
   ax = fig.add_subplot(111)

   l = 2 # grossura da linha do plot do gráfico
   fs_labels = 20 # tamanho da fonte dos rótulos dos eixos x e y.
   
   num_onda_modif = num_onda[perfil > lim_perfil]
#   s_modif = s[:, perfil > lim_perfil]
#   perfil_modif = perfil[perfil > lim_perfil]
   
#   print num_onda.shape
#   print num_onda_modif.shape
#   
#   print s.shape
#   print s_modif.shape
#   
#   print perfil.shape
#   print perfil_modif.shape
   
   for i in range(0, nl):
      ax.plot(num_onda[0], s[i,:], linewidth=l, linestyle='--')   
   
#   ax.plot(num_onda[0], s[0,:], linewidth=l, linestyle='--')
#   ax.plot(num_onda[0], s[1,:], linewidth=l, linestyle='--')
#   ax.plot(num_onda[0], s[2,:], linewidth=l, linestyle='--')
#   ax.plot(num_onda[0], s[3,:], linewidth=l, linestyle='--')
#   ax.plot(num_onda[0], s[4,:], linewidth=l, linestyle='--')
#   ax.plot(num_onda[0], s[5,:], linewidth=l, linestyle='--')
#   ax.plot(num_onda[0], s[6,:], linewidth=l, linestyle='--')
   ax.plot(num_onda[0], perfil[0], linewidth=l, linestyle='-')
   ax.plot([num_onda_cm, num_onda_cm], [0,1], linewidth=l, linestyle='--', label=u"Centro da Transição")
#   ax.plot(num_onda_cm, 0.5, "o", markersize=10, label=u"Centro da Transição")
#   ax.plot(num_onda_modif, s_modif[0,:], linewidth=l, label=u"pico 0")
#   ax.plot(num_onda_modif, s_modif[1,:], linewidth=l, linestyle='-', label=u"pico 1")
#   ax.plot(num_onda_modif, s_modif[2,:], linewidth=l, label=u"pico 2")
#   ax.plot(num_onda_modif, s_modif[3,:], linewidth=l, linestyle='-', label=u"pico 3")
#   ax.plot(num_onda_modif, s_modif[4,:], linewidth=l, linestyle='-', label=u"pico 4")
#   ax.plot(num_onda_modif, perfil_modif, linewidth=l, linestyle='-', label=u"soma dos picos")
   
   
   plt.suptitle(u"Simulação Desvios Isotópicos " + elem + " [$\lambda$ = " + str(lamb) + " nm ]", fontsize=20) # título do gráfico
   ax.legend(loc='upper left')
   #plt.xlabel(r"\textbf{n\'umero de onda (cm^{-1})}") # notacao latex
   ax.set_xlabel(u"número de onda [cm$^{-1}$]", fontsize=fs_labels)
   #plt.ylabel(r"\textbf{Intensidade (U. A.)}") # notacao latex
   ax.set_ylabel(u"Intensidade normalizada", fontsize=fs_labels)
   
   ax.set_ylim((0,1.2))
   ax.set_xlim((num_onda_modif[0], num_onda_modif[-1]))
   
      
   for i in range(len(picosx)):
      
      if M[0][i] == 163:
         ax.annotate(str(int(M[0][i])) + "-" + str((M_from[0][i])) + "-" + str((M_to[0][i])), xy=(picosx[i], picosy[i]), xytext=(picosx[i], 0.5), ha="center", rotation=90, fontsize=14, weight='bold')
             
      elif M[0][i] == 161:
         ax.annotate(str(int(M[0][i])) + "-" + str((M_from[0][i])) + "-" + str((M_to[0][i])), xy=(picosx[i], picosy[i]), xytext=(picosx[i], 0.8), ha="center", rotation=90, fontsize=14, weight='bold')
       
      else:
         ax.annotate(str(int(M[0][i])) + "-" + str((M_from[0][i])) + "-" + str((M_to[0][i])), xy=(picosx[i], picosy[i]), xytext=(picosx[i], picosy[i]+0.1), ha="center", rotation=90, fontsize=14, weight='bold')
#      ax.annotate(str(int(M[0][i])) + "-" + str((M_from[0][i])) + "-" + str((M_to[0][i])), xy=(picosx[i], picosy2[i]), xytext=(picosx[i], picosy2[i]+0.02), ha="center", rotation=90, fontsize=14, weight='bold')
   
#   ax.annotate("C.T.", xy=(num_onda_cm, 1), xytext=(num_onda_cm, 1), fontsize=14);
   
   
   # Essas duas linhas de código são respnsáveis por mudar a maneira como os números aparecem 
   # no eixo x, fazendo com que tenhamos floats mostrados com 3 casas decimais.
   locs,labels = plt.xticks()
   plt.xticks(locs, map(lambda x: "%.3f" % x, locs))
   
   
   plt.savefig(nome_fig + ".png") # salva uma figura do plot
   plt.show() # apresenta o plot na tela do computador


#---------------------------------------------------------------------
def plot_Casimir(h, filename):


   [nl, nc] = h.shape
   gama = 2 # largura da lorentziana à meia altura em MHz (FWHM = 2 * HWHM)   
   start = min(h[:,0])
   stop = max(h[:,0])
   
   
   h[:,1] = h[:,1] / max(h[:,1])
   
   if start < 0:
      start = 1.2 * start
   else:
      start = 0.8 * start
   
   if stop < 0:
      stop = 0.8 * stop
   else:
      stop = 1.2 * stop
   
   
   x = np.linspace(start, stop, num=1000)
   y = np.zeros(x.shape)
   
   for i in range(nl):
      x0 = h[i, 0]
      y0 = h[i, 1]
      y_temp = y0 * ((gama**2 / 4) / ((x - x0)**2 + gama**2 / 4))
      y = y + y_temp
      
   y = y / max(y)
   plt.rcParams['font.family'] = 'Times New Roman' # aqui a gente seta que vai usar fonte times new roman no gráfico
   fig = plt.figure()
   ax = fig.add_subplot(111)
   
   l = 2 # grossura da linha do plot do gráfico
   fs_labels = 20 # tamanho da fonte dos rótulos dos eixos x e y.
   
   ax.plot(x, y, linewidth=l, linestyle='-') 

   
   
   plt.suptitle(u"Hiperfina Casimir") # título do gráfico
   ax.legend(loc='upper left')
   #plt.xlabel(r"\textbf{n\'umero de onda (cm^{-1})}") # notacao latex
   ax.set_xlabel(u"Delta W [MHz]", fontsize=fs_labels)
   #plt.ylabel(r"\textbf{Intensidade (U. A.)}") # notacao latex
   ax.set_ylabel(u"Intensidade", fontsize=fs_labels)
   
#   ax.set_ylim((1.2*min(x),1.2*max(x)))
#   ax.set_xlim((1.2*min(y),1.2*max(y)))
   
   
   # Essas duas linhas de código são respnsáveis por mudar a maneira como os números aparecem 
   # no eixo x, fazendo com que tenhamos floats mostrados com 3 casas decimais.
   locs,labels = plt.xticks()
   plt.xticks(locs, map(lambda x: "%.3f" % x, locs))

   for i in range(0, nl):
      ax.annotate("{:.1f}".format(h[i][2]) + "-" + "{:.1f}".format(h[i][3]), xy=(h[i][0], h[i][1]), xytext=(h[i][0], h[i][1]+0.02), ha="center", rotation=90, fontsize=14, weight='bold')
   
   mat = np.vstack((x, y))
   np.savetxt(filename + ".csv", np.transpose(mat), fmt='%.18e', delimiter=',', newline='\n', header='', footer='', comments='# ')
   
   plt.savefig(filename + ".png") # salva uma figura do plot
   plt.show() # apresenta o plot na tela do computador
   
#---------------------------------------------------------------------
def setDyIni_16693():
#  valores para o Dy:
   global I, J, J_sup
   I = 5.0/2.0
   J = 8.0
   J_sup = 7.0

#---------------------------------------------------------------------
def setDyIni_16733():
#  valores para o Dy:
   global I, J, J_sup
   I = 5.0/2.0
   J = 8.0
   J_sup = 8.0

#---------------------------------------------------------------------
def setHiperfinaDy161_0_16693():
#  valores para o Isótopo 161 do Dy na transição 0 para 16693,87 cm-1:
   global Ai, Bi, Af, Bf
   Ai = -115.8     # 0 cm-1
   Bi = 1102.0    # 0 cm-1

   Af = -64.3    # 16693,87 cm-1
   Bf = 898.0     # 16693,87 cm-1
#---------------------------------------------------------------------
def setHiperfinaDy163_0_16693():
#  valores para o Isótopo 163 do Dy na transição 0 para 16693,87 cm-1:
   global filename
   filename = "hiperfina_Dy_16693.csv"
   global Ai, Bi, Af, Bf
   Ai = 162.9     # 0 cm-1
   Bi = 1150.0    # 0 cm-1

   Af = 90.0    # 16693,87 cm-1
   Bf = 947.0     # 16693,87 cm-1

#---------------------------------------------------------------------
def setHiperfinaDy161_0_16733():
#  valores para o Isótopo 161 do Dy na transição 0 para 16733,2 cm-1:
   global Ai, Bi, Af, Bf
   Ai = -115.8     # 0 cm-1
   Bi = 1102.0    # 0 cm-1

   Af = -88.6    # 16733,2 cm-1
   Bf = 1400.0     # 16733,2 cm-1

#---------------------------------------------------------------------
def setHiperfinaDy163_0_16733():
#  valores para o Isótopo 163 do Dy na transição 0 para 16733,2 cm-1:
   global filename 
   filename = "hiperfina_Dy_16733.csv"
   global Ai, Bi, Af, Bf
   Ai = 162.9     # 0 cm-1
   Bi = 1150.0    # 0 cm-1

   Af = 124.0    # 16733,2 cm-1
   Bf = 1484.0     # 16733,2 cm-1

#---------------------------------------------------------------------
def main():
   
   # adicionando o diretório corrente tanto para o scilab quanto para o octave
   sci.getd(os.getcwd())
   oc.addpath(os.getcwd())

   # Aqui nós pedimos para que os prints dos programas do scilab e do octave sejam exibidos no terminal do python
#   oc.logger.setLevel(logging.DEBUG)
   oc.logger.setLevel(logging.INFO)
#   sci.logger.setLevel(logging.INFO)
   
   # setando os valores de I, J, J_sup, Ai, Af, Bi, Bf
   setDyIni_16733()
   setHiperfinaDy161_0_16733()
   
#   setDyIni_16693()
#   setHiperfinaDy161_0_16693()
   
#   print "I = " + str(I) + ";" 
#   print "J = " + str(J) + ";"
#   print "J_sup = " + str(J_sup) + ";"
#   print "Ai = " + str(Ai) + ";"
#   print "Af = " + str(Af) + ";"
#   print "Bi = " + str(Bi) + ";"
#   print "Bf = " + str(Bf) + ";"
   
   # h161 é a tabela hiperfina do isotopo 161
   h161 = sci.callHiperfina(I, J, J_sup, Ai, Af, Bi, Bf)


# Essa linha desloca o zero pra esquerda, mas a eq. de Casimir já coloca o zero no centro de massa   
# Então não precisa dessa linha.
#   h161[:, 0] = h161[:,0] - min(h161[:,0])


   plot_Casimir(h161, "Dy_h161_16733")
#   plot_Casimir(h161, "Dy_h161_16693")
   setHiperfinaDy163_0_16733()
#   setHiperfinaDy163_0_16693()
   h163 = sci.callHiperfina(I, J, J_sup, Ai, Af, Bi, Bf)
   
   
# Essa linha desloca o zero pra esquerda, mas a eq. de Casimir já coloca o zero no centro de massa
# Então não precisa dessa linha.
#   h163[:, 0] = h163[:,0] - min(h163[:,0])

   plot_Casimir(h163, "Dy_h163_16733")
#   plot_Casimir(h163, "Dy_h163_16693")
   
   [nu, nu_cm, num_onda, num_onda_cm, perfil, s, nome_fig, M,  M_from, M_to, lamb, elem, lim_perfil] = oc.main(h161, h163)
   
   
# debug retorno do oct2py:   
#   print "##########################"
#   print "num_onda type: " + str(type(num_onda))
#   print "num_onda_cm type: " + str(type(num_onda_cm))
#   print "s type: " + str(type(s))
#   print np.where(s[0,:] == s[0,:].max())
#   print num_onda[0][np.where(s[0,:] == s[0,:].max())]
#   print perfil.shape
#   print perfil.size
#   print perfil[0,1]
#   print "##########################"
   
   plot(nu, nu_cm, num_onda, num_onda_cm, perfil, s, nome_fig, M,  M_from, M_to, lamb, elem, lim_perfil)

   mat = np.vstack((nu, num_onda, perfil, s))
   
#   print mat.shape 
#   print len(nu_cm)
#   print len(mat)
# Para salvar as matrizes para arquivo terei de usar esse comando a seguir.
# documentação em https://docs.scipy.org/doc/numpy/reference/generated/numpy.savetxt.html
   np.savetxt(filename, np.transpose(mat), fmt='%.18e', delimiter=',', newline='\n', header='', footer='', comments='# ')
#---------------------------------------------------------------------
if __name__ == "__main__":
   main()
