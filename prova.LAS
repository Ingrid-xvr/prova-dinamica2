//Prova Semestral

//Portico de 3 pavimentos com massa concentrada em cada pavimento e com 5m na horizontal(adotado, 
//uma vez que que os DOF's influenciam apenas os pilares) e 3m na vertical (dado da questao)

//Nos do portico
Int XY = [ 0	0    
	     5	0             
           0	3    
           5	3    
           0	6    
	     5	6   
	     0	9    
	     5	9 ]    
            
//Elementos do portico para serem unidos depois
Int Beam_EL = [1	3 //pilar esquerda 1 pav       
		   2	4 //pilar direita 1 pav
 		   3	4 //viga 1 pav
		   3	5 //pilar esquerda 2 pav
 		   4	6 //pilar direita 2 pav
		   5	6 //viga 2 pav
		   5	7 //pilar esquerda 3 pav
		   6	8 //pilar direita 3 pav
		   7	8] //viga 3 pav
		   
//DOF's - translacional na horizontal (x,y,rotacao)
Int EQ = [ 0 0 0    //no 1
           0 0 0    //no 2          
           1 0 0    //no 3
           1 0 0    //no 4
           2 0 0    //no 5
           2 0 0    //no 6
           3 0 0    //no 7
           3 0 0 ]  //no 8
                  			 
//Locacao das vigas na matriz
Int Beam_LM = Beam_Make_LM(Beam_EL,EQ)

// Propriedades [A As I E Nu Dn] k/2 = 12EI/l^3, portanto, E = kl^3/24I = (10^7)(3)^3/24(1)=11250000
Dbl PC = [1.0 0.0 1.0 11250000 0.0 0.0]
Dbl PB = [1.0 0.0 1.0 1.0 0.0 0.0]

// Beams stiffness matrices
Dbl KC1 = Beam_Make_K(1,XY,Beam_EL,PC,"Global")
Dbl KC2 = Beam_Make_K(2,XY,Beam_EL,PC,"Global")
Dbl KB1 = Beam_Make_K(3,XY,Beam_EL,PB,"Global")
Dbl KC3 = Beam_Make_K(4,XY,Beam_EL,PC,"Global")
Dbl KC4 = Beam_Make_K(5,XY,Beam_EL,PC,"Global")
Dbl KB2 = Beam_Make_K(6,XY,Beam_EL,PB,"Global")
Dbl KC5 = Beam_Make_K(7,XY,Beam_EL,PC,"Global")
Dbl KC6 = Beam_Make_K(8,XY,Beam_EL,PC,"Global")
Dbl KB3 = Beam_Make_K(9,XY,Beam_EL,PB,"Global")

// Matriz de Rigidez Global 3x3 (3 DOF's)
Dbl K = Zero(3)
Assemble(K,KC1,Beam_LM,1) //elemento 1
Assemble(K,KC2,Beam_LM,2) //elemento 2
Assemble(K,KB1,Beam_LM,3) //elemento 3
Assemble(K,KC3,Beam_LM,4) //elemento 4 
Assemble(K,KC4,Beam_LM,5) //elemento 5
Assemble(K,KB2,Beam_LM,6) //elemento 6 
Assemble(K,KC5,Beam_LM,7) //elemento 7
Assemble(K,KC6,Beam_LM,8) //elemento 8
Assemble(K,KB3,Beam_LM,9) //elemento 9

// Matriz de Massa - concentrada nos pavimentos (20000 kg) 
Dbl M = [20000	0		0		
		0	20000		0		
		0	0		20000]		
		
// Analise modal
Dbl {Phi,W2} = Jacobi(M,K)
Dbl W = Sqrt(W2)
Dbl Periods= W@^(-1) @* 2*3.141592654

Dbl Mgen = t(Phi)*M*Phi //Massa modal (ou massa generalizada para normalizar Phi)
Dbl Kgen = t(Phi)*K*Phi //Rigidez modal (rigidez generalizada)
Dbl Phi_T= t(Phi)

//Dados do enunciado
Dbl u0=[0.01
	  0.02
	  0.03]

Dbl z0 = Phi_T*M*u0

//Sem amortecimento 
Dbl pi = Pi
Dbl tempo = 2*pi/W(1,1)
Dbl aux = W*tempo
Dbl cos = Cos(aux)

Dbl zt1 = Col(z0,1)@*Col(cos,1)
Dbl Ut1 = Phi*zt1

//Amortecimento de 2% - Valores calculados no Excel
Dbl zt2 = [4.640497341
	    -0.055912491
	     0.039117926]

Dbl Ut2 = Phi*zt2

