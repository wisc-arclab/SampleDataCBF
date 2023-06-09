*  NLP written by GAMS Convert at 07/19/01 13:40:20
*  
*  Equation counts
*     Total       E       G       L       N       X
*        16      15       0       1       0       0
*  
*  Variable counts
*                 x       b       i     s1s     s2s      sc      si
*     Total    cont  binary integer    sos1    sos2   scont    sint
*        17      17       0       0       0       0       0       0
*  FX     0       0       0       0       0       0       0       0
*  
*  Nonzero counts
*     Total   const      NL     DLL
*        45      33      12       0
*
*  Solve m using NLP minimizing objvar;


Variables  x1,x2,objvar,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17;

Positive Variables x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17;

Equations  e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16;


e1..  - 3*x1 - 3*x2 - objvar + 2*x4 + 2*x5 =E= 60;
* e1.. x1 + objvar =E= 0;
  
e2..    x1 - 2*x2 + x4 + x5 =L= 40;

e3..    2*x1 - x4 + x6 =E= -10;

e4..    2*x2 - x5 + x7 =E= -10;

e5..  - x1 + x8 =E= 10;

e6..    x1 + x9 =E= 20;

e7..  - x2 + x10 =E= 10;

e8..    x2 + x11 =E= 20;

e9.. x6*x12 =E= 0;

e10.. x7*x13 =E= 0;

e11.. x8*x14 =E= 0;

e12.. x9*x15 =E= 0;

e13.. x10*x16 =E= 0;

e14.. x11*x17 =E= 0;

e15..    2*x1 - 2*x4 + 2*x12 - x14 + x15 =E= -40;

e16..    2*x2 - 2*x5 + 2*x13 - x16 + x17 =E= -40;

* set non default bounds

*
* Lower and upper bounds of variables have been computed by SparsePOP
*
x1.lo = -10;
x2.up = 5;
x2.lo = -10;
x2.up = 5;

x4.up = 22.5; 
x5.up = 25; 
x6.up = 20; 
x7.up = 20; 
x8.up = 12.5; 
x9.up = 30; 
x10.up = 15; 
x11.up = 30; 
x12.up = 0.75; 
x13.up = 0.2; 
x14.up = 20; 
x15.up = 0.3; 
x16.up = 20; 
x17.up = 0.1; 

* set non default levels

x1.l = -8; 
x2.l = -8; 
x4.l = 1; 

* set non default marginals


Model m / all /;

m.limrow=0; m.limcol=0;

$if NOT '%gams.u1%' == '' $include '%gams.u1%'

Solve m using NLP minimizing objvar;