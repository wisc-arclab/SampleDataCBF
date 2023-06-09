*  NLP written by GAMS Convert at 08/29/02 12:53:04
*  
*  Equation counts
*     Total       E       G       L       N       X       C
*         2       2       0       0       0       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*         3       3       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*         5       3       2       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,objvar;

Positive Variables x1,x2;

Equations  e1,e2;


e1..  -x2^2 + 7*x2 + 12*x1 + objvar =E= 0;

e2..  - 2*x1^4 - x2 =E= -2;

* set non default bounds

x1.up = 2; 
x2.up = 3; 

* set non default levels

x1.l = 0.7175; 
x2.l = 1.47; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;
