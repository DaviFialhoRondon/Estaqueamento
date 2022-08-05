%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  ENTRADA  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Número de estacas
numEstacas = 16;

#Coordenadas topo estacas (m)
#              1        2      3        4      5       6      7       8      9      10      11      12     13      14     15      16
 corXEst = [ -3.00 ; -3.00 ; -1.50 ; -1.50 ; 0.00 ;  0.00 ; 1.50 ;  1.50 ; 3.00 ;  3.00 ; -1.50 ; -1.50 ; 0.00 ;  0.00 ; 1.50 ;  1.50 ] ;
 corYEst = [  0.75 ; -0.75 ;  2.25 ; -2.25 ; 2.25 ; -2.25 ; 2.25 ; -2.25 ; 0.75 ; -0.75 ;  0.75 ; -0.75 ; 0.75 ; -0.75 ; 0.75 ; -0.75 ] ;

#Inclinação vertical estacas (deg)
#         1    2    3    4    5    6    7    8   9    10   11  12  13  14  15  16
beta = [ 10 ; 10 ; 10 ; 10 ; 10 ; 10 ; 10 ; 10 ; 10 ; 10 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ] ;

#Angulo horizontal estacas (deg)
#         1     2     3     4     5    6    7   8   9  10   11  12  13  14  15  16
gama = [ 180 ; 180 ; 180 ; 180 ; 90 ; 270 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ] ;

#Cargas na horigem do bloco
#          FX(kN)  FY(kN)  FZ(kN)  MX(kNm)    MY(kNm)   MZ(kNm)
cargas = [  185  ;  55  ;  -8350 ;   240   ;   1620   ;    0    ] ;


printf(" \t\t      X  \t  Y \t Beta \tGama \n");
for est = 1 : numEstacas
     printf(" estaca %2.0i  -  %10.2f  %10.2f  %5.0f   %5.0f\n",est,corXEst(est),corYEst(est),beta(est),gama(est));

endfor



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  MATRIZ ROTAÇÃO E TRANSLAÇÃO  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mrot] = Rotation(B,G)

     Mrot(1,1) =  sind(B) * cosd(G);
     Mrot(1,2) =  sind(B) * sind(G);
     Mrot(1,3) = -cosd(B)          ;

     Mrot(2,1) =  cosd(B) * cosd(G);
     Mrot(2,2) =  cosd(B) * sind(G);
     Mrot(1,3) =  sind(B)          ;

     Mrot(3,1) =  sind(G) ;
     Mrot(3,2) = -cosd(G) ;
     Mrot(3,3) =  0       ;

endfunction

function [Mtra] = Translation(X,Y)

     Mtra(1,1) =  1 ;
     Mtra(1,2) =  0 ;
     Mtra(1,3) =  0 ;
     Mtra(1,4) =  0 ;
     Mtra(1,5) =  0 ;
     Mtra(1,6) = -Y ;

     Mtra(2,1) =  0 ;
     Mtra(2,2) =  1 ;
     Mtra(2,3) =  0 ;
     Mtra(2,4) =  0 ;
     Mtra(2,5) =  0 ;
     Mtra(2,6) =  X ;

     Mtra(3,1) =  0 ;
     Mtra(3,2) =  0 ;
     Mtra(3,3) =  1 ;
     Mtra(3,4) =  Y ;
     Mtra(3,5) = -X ;
     Mtra(3,6) =  0 ;

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  SISTEMA  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kLocal  = zeros(3,3) ; #Matriz de rigidez local
kGlobal = zeros(6,6) ; #Matriz de rotação global

matRotation    = zeros(3,3) ; #Matriz de rotação
matTranslation = zeros(3,6) ; #Matriz translação

#Matriz de rigidez dos elementos (EA/L)
kLocal(1,1) = 1;

#Formação da matriz rigidez global
for estaca = 1 : numEstacas

     matRotation    = Rotation    ( beta(estaca)    , gama(estaca)    );
     matTranslation = Translation ( corXEst(estaca) , corYEst(estaca) );

     kGlobal = kGlobal + matTranslation' * matRotation' * kLocal * matRotation * matTranslation ;

endfor

#Resolução do sistema K * U = F
     uGlobal = kGlobal^(-1) * cargas ; #deslocamento global

uGlobal

#Cálculo dos esforços em cada estaca
for estaca = 1 : numEstacas

     matRotation    = Rotation    ( beta(estaca)    , gama(estaca)    );
     matTranslation = Translation ( corXEst(estaca) , corYEst(estaca) );

     uLocal = matRotation * matTranslation * uGlobal  ;

     printf(' estaca[%2i] = %5.3f kN\n', estaca , uLocal(1));
endfor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  SAÍDAS  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








