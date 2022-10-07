#****************************************************
#                                                   *
#       MODELO DE ISING EM DUAS DIMENSOES           *
#                                                   *
#****************************************************

# Programa para calculo das caracteristicas fisicas de um sistema de spins bidimensional
# descrito pelo Modelo de Ising usando Monte Carlo e o algoritmo de Metropolis.

tempo1 = time(); # para calcular o tempo que o programa demora para rodar

# Variaveis a serem usadas dentros das funcoes:
global n_spins T H J n_var indice_H_nulo



#***********PARTE QUE PODE SER ALTERADA: ************************************************

# Parametros iniciais:
n_spins = [40];                           # numero de spins em cada dimensão, vetor
T = [1:0.2:2 2.1:0.05:2.3 2.4:0.2:3.5];   # temperaturas, vetor
H = [-1 1];                                  # valores de campo magnetico, vetor
n_var = [2000];                           # numero de varreduras, vetor
J = 1;                                    # unidade de energia, escalar

# Zoom ao redor de T=2.2:
temperatura_ini = 2.1;  # temperatura inicial que voce quer dar zoom
temperatura_fim = 2.3;  # temperatura final que voce quer dar zoom

#***********FIM**************************************************************************





#********************TEMPO DO PROGRAMA: *************************************************

# Duracao aproximada do programa:
# Caso veja que va demorar muito tempo pode parar o programa e escolher um
# numero menor de passos em cada elemento. 
# Lembrando que (em segundos):
# tempo de execução =~ 0.00018151 * sum(n_spins.^2) * n_var * length(H) * length(T)  
segundos = 0.00018151 * sum(n_spins.^2) * sum(n_var) * length(H) * length(T);
disp ('Tempo aproximado de execucao:')
disp('(hh:mm:ss)')
disp(datestr(segundos/(24*60*60), 'HH:MM:SS'))

#********************AJUSTES INICIAIS: **************************************************


# Arrays para os calculos:
spins = ones(n_spins(1));                                    # sitio de spins, matriz
m = zeros(length(T),length(H),length(n_spins),length(n_var));# magnetizacao media por sitio a cada valor de H,de T e de n_spins, tensor
E = zeros(length(T),1);                                      # energia media por sitio a cada T, vetor
C = zeros(length(T),1);                                      # calor especifico a cada T, vetor

# E e C sao calculados sempre para o mesmo H num dado tamanho de grid e com determinado n_var

# Ajuste do indice do vetor H para H=0:
indice_H_nulo = 1; # se não tiver H=0 usa o primeiro valor
for b=1:length(H)
  if H(b)== 0
    indice_H_nulo = b;
  endif
endfor

# Ajuste dos indices do vetor T para as temperaturas iniciais e finais escolhidas:
indice_t_ini = indice_t_fim = 1;
for a=1:length(T)
  if abs( temperatura_ini-T(a) ) < abs( T(indice_t_ini)-temperatura_ini ) 
    indice_t_ini = a;
  endif
  if abs( temperatura_fim-T(a) ) < abs( T(indice_t_fim)-temperatura_fim )  
    indice_t_fim = a;
  endif
endfor


#********************FUNCOES: ***********************************************************

# FUNCAO 1
# Funcao que calcula a energia total do sistema:
# recebe o grid de spins e o valor do campo H
function f = energia(spins,h)
  
  global J
  
  # largura do grid de spins:
  n_spins = length(spins(:,1));
  
  # energia no canto inferior direito do grid:
  energia = -spins(n_spins,n_spins)*(J*(spins(n_spins,1)+spins(1,n_spins)) + h);
  
  # energia na linha inferior do grid:
  ii=n_spins;
  for j=1:n_spins-1
    energia = energia-spins(ii,j)*(J*(spins(ii,j+1)+spins(1,j)) + h);
  endfor
  
  # energia na lateral direita do grid:
  j=n_spins;
  for ii=1:n_spins-1
    energia = energia-spins(ii,j)*(J*(spins(ii,1)+spins(ii+1,j)) + h);
  endfor
  
  # energia no 'meio' dao grid:
  for ii=1:n_spins-1
    for j=1:n_spins-1
      energia = energia-spins(ii,j)*(J*(spins(ii,j+1)+spins(ii+1,j)) + h);
    endfor
  endfor
  
  # retorna a energia total do sistema:
  f = energia;
  
endfunction

# FUNCAO 2
# Funcao que calcula E_flip:  
# recebe o grid de spins, a posicao do spin (i,j), o valor do campo H
function f = flip(spins,ii,j,h) 
  
  global J
  
  # largura do grid de spins:
  n_spins = length(spins(:,1));
  
  #arrumar as bordas da matriz de spins devido as condicoes periodicas de contorno:
  if ii == 1
    i_ant = n_spins;
    i_pos = ii+1;
  elseif ii == n_spins
    i_ant = ii-1;
    i_pos = 1;
  else
    i_ant = ii-1;
    i_pos = ii+1;
  endif
  
  if j == 1
    j_ant = n_spins;
    j_pos = j+1;
  elseif j == n_spins
    j_ant = j-1;
    j_pos = 1;
  else
    j_ant = j-1;
    j_pos = j+1;
  endif
  
  # retorna a energia E_flip:
  f = 2*spins(ii,j)*( J*(spins(i_ant,j)+spins(i_pos,j)+spins(ii,j_ant)+spins(ii,j_pos)) + 2*h );
  
endfunction

# FUNCAO 3
# Funcao que calcula m, E e C para uma dada configuracao:
# recebe o grid de spins, indices do vetor n_var, T e H  
function [mag, Ene, Cap] = calcula_config(spins,n,v,t,h)  
  
  global n_spins T H n_var indice_H_nulo
    
  m_a = zeros(n_var(v),1);      # magnetizacao media por sitio a cada varredura, vetor
  E_a = zeros(n_var(v),1);      # energia media por sitio a cada varredura, vetor  

  for varredura=1:n_var(v)

    # Algoritmo de MC - Metropolis:  
    for ii=1:n_spins(n)      # spins da cadeia, linha
      for j=1:n_spins(n)     # spins da cadeia, coluna
        
        E_flip = flip(spins,ii,j,H(h));       #calcula a energia
        if E_flip < 0
          spins(ii,j) = -spins(ii,j);         #flipa o spin
        else
          P_flip = exp(-E_flip/T(t));         #probabilidade
          r = rand();                         #numero sorteado
          if r < P_flip 
            spins(ii,j) = -spins(ii,j);       #flipa o spin
          endif
        endif
        
      endfor                 # acabam colunas
    endfor                   # acabam linhas  

    # calcula m_alpha (magnetizacao) apos cada varredura:
    m_a(varredura) = sum(sum(spins)) / n_spins(n)^2;  
    if h==indice_H_nulo # indice do vetor H no qual H=0 (se H tiver esse valor, se nao =H(1) )
      # energia media por spin do sistema apos cada varredura:
      E_a(varredura) = energia(spins,H(h)) / n_spins(n)^2;
    endif

  endfor                    # acabam varreduras 

  # calcula <m> :
  mag = sum(m_a) / n_var(v);
  
  # caso H!=0 nao se usam esses valores:
  Ene = Cap = 0;  
  
  # calcula E e C para H=0 (se H tiver esse valor, se nao para H(1) )
  if h==indice_H_nulo 
    # calcular a energia media por sitio:
    Ene = sum(E_a) / n_var(v);
    # calcula o calor especifico:
    Cap = var(E_a) / (T(t)^2);
  endif
  
endfunction  
 
 
#********************PRINCIPAL: *********************************************************
  
# Calculo do arranjo dos spins e magnetizacao:
for n=1:length(n_spins)        # valores do tamanho do grid, n = indices do vetor n_spins
  for v=1:length(n_var)        # varreduras, v = indices do vetor n_var
    for t=1:length(T)          # valores de temperatura, t = indices do vetor T
      spins = ones(n_spins(n),n_spins(n)); # gera um sitio de spins no novo tamanho para cada T    
      for h=1:length(H)        # valores de campo magnetico, b = indices do vetor H
        
        if h==indice_H_nulo
          [mag, Ene, Cap] = calcula_config(spins,n,v,t,h); # recebe valores de m, E e C para uma dada configuracao
          m(t,h,n,v) = mag;
          E(t) = Ene;
          C(t) = Cap;
        else
          [mag, Ene, Cap] = calcula_config(spins,n,v,t,h); # recebe valores de m
          m(t,h,n,v) = mag;
        endif
        
      endfor                   # acabam valores de campo magnetico  
    endfor                     # acabam valores de temperatura
  endfor                       # acabam valores de numero de varreduras
endfor                         # acabam valores de tamanho do grid


#********************SALVAR OS DADOS: ***************************************************

# Salvar as variaveis para poderem ser lidas posteriormente:
tempo = ctime (time ());
for r=[4 8 11 14 17 20]
  tempo(r) = "_";
endfor
arquivo = strcat(tempo(9:11),tempo(5:8),tempo(12:19),'.txt');
save(arquivo,'n_var','n_spins','T','H','E','C','J','m','indice_t_ini','indice_t_fim','indice_H_nulo');
# o arquivo e salvo com a data e a hora em que foi feito.


#********************GRAFICOS: ********************************************************** 

# ultimos valores com que E e C foram calculados:  
indice_spins = length(n_spins);
indice_var = length(n_var);   
  
# GRAFICO 1  
# Plot da Magnetizacao media por sitio em funcao de T com H=0 (se H tiver esse valor, se nao =H(1) )
# m x T para um H (1 linha no grafico)
figure;
plot(T,abs(m(:,indice_H_nulo,1,1)),'o-b');                     # m para todas as T, para o primeiro valor do vetor n_spins e H=0 (se H tiver esse valor, se nao =H(1) )
titulo = strcat('Magnetizacao media x Temperatura. H=',num2str(H(indice_H_nulo))); # gera um texto para o titulo
title(titulo,'FontSize',20);                                   # coloca um titulo
xlabel('Temperatura','FontSize',20);                           # coloca uma legenda no eixo x
ylabel('<m> por sitio','FontSize',20);                         # coloca uma legenda no eixo y
legenda = strcat("H=",num2str(H(indice_H_nulo)));              # gera um texto pra legenda
legenda2 = legend(legenda,'location','northeast');             # adiciona uma legenda
set(legenda2,'FontSize',20)                                    # tamanho da legenda
set(gca,'FontSize',20)                                         # tamanho dos numeros nos eixos
yl = ylim();                                                   # devolve os limites no eixo y
xl = xlim();                                                   # devolve os limites no eixo x
texto = strcat('Grid = ',num2str(n_spins(1)),'x',num2str(n_spins(1)));# gera uma string
texto2 = strcat('Nvar = ',num2str(n_var(1)));                  # gera uma string
text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,texto,'FontSize',16); # coloca um texto para o tamanha da cadeia de spins
text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.75,texto2,'FontSize',16);# coloca um texto para o numero de varreduras
  
  
# GRAFICO 2  
# Plot da Energia por sitio em funcao de T com H=0 (se H tiver esse valor, se nao =H(1) )
# E x T para um H (1 linha no grafico)
figure;
plot(T,E,'o-b');                                            # E para todas as T com H=0 (se H tiver esse valor, se nao =H(1) )
titulo = strcat('Energia x Temperatura. H=',num2str(H(indice_H_nulo))); # gera um texto para o titulo
title(titulo,'FontSize',20);                                   # coloca um titulo
xlabel('Temperatura','FontSize',20);                           # coloca uma legenda no eixo x
ylabel('<E> por sitio','FontSize',20);                         # coloca uma legenda no eixo y
legenda = strcat("H=",num2str(H(indice_H_nulo)));              # gera um texto pra legenda
legenda2 = legend(legenda,'location','northwest');             # adiciona uma legenda
set(legenda2,'FontSize',20)                                    # tamanho da legenda
set(gca,'FontSize',20)                                         # tamanho dos numeros nos eixos
yl = ylim();                                                   # devolve os limites no eixo y
xl = xlim();                                                   # devolve os limites no eixo x
texto = strcat('Grid = ',num2str(n_spins(indice_spins)),'x',num2str(n_spins(indice_spins)));# gera uma string
texto2 = strcat('Nvar = ',num2str(n_var(indice_var)));                  # gera uma string
text(xl(1)+(xl(2)-xl(1))*0.04,yl(1)+(yl(2)-yl(1))*0.75,texto,'FontSize',16); # coloca um texto para o tamanha da cadeia de spins
text(xl(1)+(xl(2)-xl(1))*0.04,yl(1)+(yl(2)-yl(1))*0.65,texto2,'FontSize',16);# coloca um texto para o numero de varreduras
  
  
# GRAFICO 3 --- determinar numero de varreduras n_var 
# Plot da Magnetizacao media por sitio em funcao de T com H=0 (se H tiver esse valor, se nao =H(1) )
# para varios tamanhos de n_var
# varias T e um H, varios n_var  ( length(n_var) linhas no grafico ) 
estilo = {'o-r','s-b','x-g','^-k','p-m','+-c'}; 
if length(n_var)>1
  figure;
  for n=1:length(n_var) # valores de n_var
    plot(T, abs(m(:,indice_H_nulo,1,n)), estilo{ mod(n,length(estilo))+1 } ); hold on; # plota uma linha para cada n_var, para H=0 (se H tiver esse valor, se nao =H(1) ) 
    legendCell{n} = num2str(n_var(n),'Nvar=%-d');                # cria uma legenda
  endfor
  legenda2 = legend(legendCell,'location','northeast');          # posiciona a legenda
  set(legenda2,'FontSize',20)                                    # tamanho da legenda
  titulo = strcat('Magnetizacao media x Temperatura. H=',num2str(H(indice_H_nulo))); # gera um texto para o titulo
  title(titulo,'FontSize',20);                                   # coloca um titulo
  xlabel('Temperatura','FontSize',20);                           # coloca uma legenda no eixo x
  ylabel('<m> por sitio','FontSize',20);                         # coloca uma legenda no eixo y
  set(gca,'FontSize',20)                                         # tamanho dos numeros nos eixos
  yl = ylim();                                                   # devolve os limites no eixo y
  xl = xlim();                                                   # devolve os limites no eixo x
  texto = strcat('Grid = ',num2str(n_spins(1)),'x',num2str(n_spins(1)));     # gera uma string
  text(xl(1)+(xl(2)-xl(1))*0.05,yl(1)+(yl(2)-yl(1))*0.55,texto,'FontSize',16);# coloca um texto para o numero de varreduras
endif


# GRAFICO 4 --- determinar tamanho de grid
# Plot da Magnetizacao media por sitio em funcao de T com H=0 (se H tiver esse valor, se nao =H(1) )
# para varios tamanhos de grid
# m x T para um H. isso para varios n_spins ( length(n_spins) linhas no grafico )
if length(n_spins)>1
  figure;
  for n=1:length(n_spins) # valores de n_spins
    plot(T,abs(m(:,indice_H_nulo,n,1)), estilo{ mod(n,length(estilo))+1 } ); hold on; # plota uma linha para cada n_spins, para H=0 (se H tiver esse valor, se nao =H(1) ) 
    legendCell{n} = num2str(n_spins(n),'L=%-d');                 # cria uma legenda
  endfor
  legenda2 = legend(legendCell,'location','northeast');          # posiciona a legenda
  set(legenda2,'FontSize',20)                                    # tamanho da legenda
  titulo = strcat('Magnetizacao media x Temperatura. H=',num2str(H(indice_H_nulo))); # gera um texto para o titulo
  title(titulo,'FontSize',20);                                   # coloca um titulo
  xlabel('Temperatura','FontSize',20);                           # coloca uma legenda no eixo x
  ylabel('<m> por sitio','FontSize',20);                         # coloca uma legenda no eixo y
  set(gca,'FontSize',20)                                         # tamanho dos numeros nos eixos
  yl = ylim();                                                   # devolve os limites no eixo y
  xl = xlim();                                                   # devolve os limites no eixo x
  texto2 = strcat('Nvar = ',num2str(n_var(1)));                  # gera uma string
  text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.55,texto2,'FontSize',16);# coloca um texto para o numero de varreduras
endif


# A PARTIR DAQUI JA FOI DEFINIDO n_spins E n_var !!!  

  
# GRAFICO 5 --- semelhante ao grafico 1 mas com zoom
# Plot da Magnetizacao media por sitio em funcao de T com H=0 (se H tiver esse valor, se nao =H(1) )
# Transicao de fase: zoom ao redor de T=2.2
# m x T para um H (1 linha no grafico)
figure;
# zoom:
x = real(log(T(indice_t_ini:indice_t_fim)-2.269185));
y = real(log(abs( m(indice_t_ini:indice_t_fim,indice_H_nulo,1,1) )));
# ajuste linear:
F = [log(ones(length(x), 1)), x(:)];
[p, e_var, r, p_var, fit_var] = LinearRegression (F, y);         # regressao linear
parametros_linear_e_angular_e_respectivas_incertezas = real([p, sqrt(p_var)]) # valores obtidos no ajuste
yFit = F * p;                                                    # ajuste
plot(x, y, '+b', x, yFit, '-g',...                               # plot dos dados e ajuste
     x, yFit + 1.96 * sqrt (fit_var), '--k',...
     x, yFit - 1.96 * sqrt (fit_var), '--k');
grid on;                                                         # coloca um grid no grafico
legenda = legend('dados','ajuste','+/- 95% valores do ajuste','location','northwest');
set(legenda,'FontSize',20)                                       # tamanho da legenda
titulo = strcat('Magnetizacao media x Temperatura. log-log. H=',num2str(H(indice_H_nulo))); # gera um texto para o titulo
title(titulo,'FontSize',20);                                     # coloca um titulo
xlabel('log(T-Tc)','FontSize',20);                               # coloca uma legenda no eixo x
ylabel('log(<m> por sitio)','FontSize',20);                      # coloca uma legenda no eixo y
set(gca,'FontSize',20)                                           # tamanho dos numeros nos eixos
yl = ylim();                                                     # devolve os limites no eixo y
xl = xlim();                                                     # devolve os limites no eixo x
texto = strcat('Grid = ',num2str(n_spins(1)),'x',num2str(n_spins(1))); # gera uma string
texto2 = strcat('Nvar = ',num2str(n_var(1)));                    # gera uma string
text(xl(1)+(xl(2)-xl(1))*0.75,yl(1)+(yl(2)-yl(1))*0.2,texto,'FontSize',16); # coloca um texto para o tamanha da cadeia de spins
text(xl(1)+(xl(2)-xl(1))*0.75,yl(1)+(yl(2)-yl(1))*0.1,texto2,'FontSize',16);# coloca um texto para o numero de varreduras
   
  
# GRAFICO 6  
# Plot do calor especifico em funcao de T com H=0 (se H tiver esse valor, se nao =H(1) )
# C x T para um H (1 linha no grafico)
figure;
plot(T,C,'o-b');                                            # C para todas as T, para H = 0 (se H tiver esse valor, se nao =H(1) ) 
titulo = strcat('Calor especifico x Temperatura. H=',num2str(H(indice_H_nulo))); # gera um texto para o titulo
title(titulo,'FontSize',20);                                   # coloca um titulo
xlabel('Temperatura','FontSize',20);                           # coloca uma legenda no eixo x
ylabel('Calor especifico','FontSize',20);                      # coloca uma legenda no eixo y
legenda = strcat('H=',num2str(H(indice_H_nulo)));              # gera um texto para a legenda
legenda2 = legend(legenda,'location','northeast');             # adiciona uma legenda
set(legenda2,'FontSize',20)                                    # tamanho da legenda
set(gca,'FontSize',20)                                         # tamanho dos numeros nos eixos
yl = ylim();                                                   # devolve os limites no eixo y
xl = xlim();                                                   # devolve os limites no eixo x
texto = strcat('Grid = ',num2str(n_spins(indice_spins)),'x',num2str(n_spins(indice_spins)));# gera uma string
texto2 = strcat('Nvar = ',num2str(n_var(indice_var)));                  # gera uma string
text(xl(1)+(xl(2)-xl(1))*0.05,yl(1)+(yl(2)-yl(1))*0.85,texto,'FontSize',16); # coloca um texto para o tamanha da cadeia de spins
text(xl(1)+(xl(2)-xl(1))*0.05,yl(1)+(yl(2)-yl(1))*0.75,texto2,'FontSize',16);# coloca um texto para o numero de varreduras
  
 
# GRAFICO 7 --- semelhante ao grafico 5 mas com zoom
# Plot do calor especifico em funcao de T com H=0 (se H tiver esse valor, se nao =H(1) )
# Transicao de fase: zoom ao redor de T=2.2
# C x T para um H (1 linha no grafico)
figure;
# zoom:
x = real(log(T(indice_t_ini:indice_t_fim)-2.269185));
y = real(log(C(indice_t_ini:indice_t_fim)));
# ajuste linear:
F = [ones(length(x), 1), x(:)];
[p, e_var, r, p_var, fit_var] = LinearRegression (F, y);         # regressao linear
parametros_linear_e_angular_e_respectivas_incertezas = real([p, sqrt(p_var)]) # valores obtidos no ajuste
yFit = F * p;                                                    # ajuste
plot(x, y, '+b', x, yFit, '-g',...                               # plot dos dados e ajuste
     x, yFit + 1.96 * sqrt (fit_var), '--k',...
     x, yFit - 1.96 * sqrt (fit_var), '--k');
grid on;                                                         # coloca um grid no grafico
legenda = legend('dados','ajuste','+/- 95% valores do ajuste','location','northeast');
set(legenda,'FontSize',20)                                       # tamanho da legenda
titulo = strcat('Calor especifico x Temperatura. log-log. H=',num2str(H(indice_H_nulo))); # gera um texto para o titulo
title(titulo,'FontSize',20);                                     # coloca um titulo
xlabel('log(T-Tc)','FontSize',20);                               # coloca uma legenda no eixo x
ylabel('log(calor especifico)','FontSize',20);                   # coloca uma legenda no eixo y
set(gca,'FontSize',20)                                           # tamanho dos numeros nos eixos
yl = ylim();                                                     # devolve os limites no eixo y
xl = xlim();                                                     # devolve os limites no eixo x
texto = strcat('Grid = ',num2str(n_spins(indice_spins)),'x',num2str(n_spins(indice_spins))); # gera uma string
texto2 = strcat('Nvar = ',num2str(n_var(indice_var)));                    # gera uma string
text(xl(1)+(xl(2)-xl(1))*0.4,yl(1)+(yl(2)-yl(1))*0.9,texto,'FontSize',16); # coloca um texto para o tamanha da cadeia de spins
text(xl(1)+(xl(2)-xl(1))*0.4,yl(1)+(yl(2)-yl(1))*0.8,texto2,'FontSize',16);# coloca um texto para o numero de varreduras
   
  
# GRAFICO 8
# Plot da Magnetizacao media por sitio em funcao de T para varios H
# m x T para varios valores de H ( length(H) linhas no grafico )
if length(H)>1
  figure;
  for c=1:length(H) # valores de H
    plot(T,m(:,c,1,1), estilo{ mod(c,length(estilo))+1 } );hold on;# plota uma linha para cada valor de H, sempre com o primeiro valor do vetor n_spins
    legendCell{c} = num2str(H(c),'N=%-d');                       # cria uma legenda
  endfor
  legenda2 = legend(legendCell,'location','northwest');          # adiciona uma legenda
  set(legenda2,'FontSize',20)                                    # tamanho da legenda
  title('Magnetizacao media x H para varios valores de T','FontSize',20);# coloca um titulo
  xlabel('H','FontSize',20);                                     # coloca uma legenda no eixo x
  ylabel('<m> por sitio','FontSize',20);                         # coloca uma legenda no eixo y
  set(gca,'FontSize',20)                                         # tamanho dos numeros nos eixos
  yl = ylim();                                                   # devolve os limites no eixo y
  xl = xlim();                                                   # devolve os limites no eixo x
  texto = strcat('Grid = ',num2str(n_spins(1)),'x',num2str(n_spins(1))); # gera uma string
  texto2 = strcat('Nvar = ',num2str(n_var(1)));                  # gera uma string
  text(xl(1)+(xl(2)-xl(1))*0.4,yl(1)+(yl(2)-yl(1))*0.75,texto,'FontSize',16); # coloca um texto para o tamanha da cadeia de spins
  text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.65,texto2,'FontSize',16);# coloca um texto para o numero de varreduras
endif  
  
  
#********************MARCA O TEMPO TOTAL DO PROGRAMA: ***********************************
  
tempo2 = time(); #tempo que o programa roda
#tempo que o programa demorou:
disp('')
disp ('Tempo que o programa demorou:')
disp('(hh:mm:ss)')
disp(datestr((tempo2-tempo1)/(24*60*60), 'HH:MM:SS'))
  
  
#********************FIM ****************************************************************
  

  
  