clc
clear all
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  ENTRADA  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#numero de tipos de estacas
numTipos = 2 ;

#Tipo para cada estaca
#               1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
tiposEstaca = [ 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 2 , 2 ] ;

#Área da seção para cada tipo (m²)
tiposArea = [ 1 , 1 ] ;

#Área da seção para cada tipo (kN/m²)
tiposElas = [ 280000 , 280000 ] ;

#Comprimento de cada tipo (m)
tiposComp = [ 12 , 12 ] ;

#Número de estacas
numEstacas = 16;

#Coordenadas topo estacas (m)
#              1        2      3        4      5       6      7       8      9      10      11      12     13      14     15      16
 corXEst = [ -3.00 ; -3.00 ; -1.50 ; -1.50 ; 0.00 ;  0.00 ; 1.50 ;  1.50 ; 3.00 ;  3.00 ; -1.50 ; -1.50 ; 0.00 ;  0.00 ; 1.50 ;  1.50 ] ;
 corYEst = [  0.75 ; -0.75 ;  2.25 ; -2.25 ; 2.25 ; -2.25 ; 2.25 ; -2.25 ; 0.75 ; -0.75 ;  0.75 ; -0.75 ; 0.75 ; -0.75 ; 0.75 ; -0.75 ] ;

#Inclinação vertical estacas (deg)
#          1      2      3      4      5      6      7      8      9     10    11     12    13    14    15    16
beta = [ 10.0 ; 10.0 ; 10.0 ; 10.0 ; 10.0 ; 10.0 ; 10.0 ; 10.0 ; 10.0 ; 10.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ] ;

#Angulo horizontal estacas (deg)
#          1       2       3       4      5      6       7     8     9     10    11    12    13    14    15    16
gama = [ 180.0 ; 180.0 ; 180.0 ; 180.0 ; 90.0 ; 270.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ; 0.0 ] ;

#Cargas na horigem do bloco
#           FX(kN)    FY(kN)    FZ(kN)    MX(kNm)      MY(kNm)     MZ(kNm)
cargas = [  185.0  ;  55.0  ;  -8350.0 ;   240.0   ;   1620.0   ;    0.0    ] ;


rigidezEstaca = zeros(numEstacas,1);

for estaca = 1 : numEstacas

     area  = tiposArea( tiposEstaca(estaca) ) ;
     elast = tiposElas( tiposEstaca(estaca) ) ;
     comp  = tiposComp( tiposEstaca(estaca) ) ;

     rigidezEstaca(estaca) = area * elast / comp ;

endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  MATRIZ ROTAÇÃO E TRANSLAÇÃO  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mrot] = Rotation(B,G)

     Mrot(1,1) =  sind(B) * cosd(G) ;
     Mrot(1,2) =  sind(B) * sind(G) ;
     Mrot(1,3) = -cosd(B)           ;
     Mrot(2,1) =  cosd(B) * cosd(G) ;
     Mrot(2,2) =  cosd(B) * sind(G) ;
     Mrot(2,3) =  sind(B)           ;
     Mrot(3,1) =  sind(G)           ;
     Mrot(3,2) = -cosd(G)           ;

endfunction

function [Mtra] = Translation(X,Y)

     Mtra(1,1) =  1.0 ;
     Mtra(1,6) = -Y   ;
     Mtra(2,2) =  1.0 ;
     Mtra(2,6) =  X   ;
     Mtra(3,3) =  1.0 ;
     Mtra(3,4) =  Y   ;
     Mtra(3,5) = -X   ;

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  ANÁLISE  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kLocal  = zeros(3,3) ; #Matriz de rigidez local
kGlobal = zeros(6,6) ; #Matriz de rotação global

matRotation    = zeros(3,3) ; #Matriz de rotação
matTranslation = zeros(3,6) ; #Matriz translação

#Matriz de rigidez dos elementos (EA/L)
kLocal         = zeros(3,3) ; #Matriz de rigidez local

#Formação da matriz rigidez global
for estaca = 1 : numEstacas

     matRotation    = Rotation    ( beta(estaca)    , gama(estaca)    ) ;
     matTranslation = Translation ( corXEst(estaca) , corYEst(estaca) ) ;

     kLocal(1,1)    = 1.0 * rigidezEstaca(estaca) ;

     kGlobal = kGlobal + matTranslation' * matRotation' * kLocal * matRotation * matTranslation ;

endfor


#Resolução do sistema K * U = F
     uGlobal = inv(kGlobal) * cargas ; #deslocamento global

#Cálculo dos esforços em cada estaca (kN)
normais = zeros(numEstacas,1) ;

for estaca = 1 : numEstacas

     matRotation    = Rotation    ( beta(estaca)    , gama(estaca)    );
     matTranslation = Translation ( corXEst(estaca) , corYEst(estaca) );

     uLocal =  matRotation * matTranslation * uGlobal ;

     normais(estaca) = rigidezEstaca(estaca) * uLocal(1)  ;
endfor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  SAÍDAS  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf('\n');
printf('\t\t\t ENTRADA DE DADOS\n\n');

printf("   Estaca        X       Y      Beta    Gama   Rigidez Axial\n");
for est = 1 : numEstacas

     printf(" estaca %2.0i  -  %5.2f | %5.2f | %4.0f  | %4.0f  |   %5.2f   |\n",est,corXEst(est),corYEst(est),beta(est),gama(est),rigidezEstaca(estaca));

endfor


printf('\n\n');


printf("   Estaca        Normal \n");
for estaca = 1 : numEstacas
     printf(' estaca %2i =   %5.3f kN\n', estaca , normais(estaca));
endfor





