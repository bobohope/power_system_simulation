
disp('9bus3Machine')
disp('------------------------------------');
global bus line mac_con pss_con exc_con


%    b#  Vm			Va			Pgen	Qgen	Pl		Ql		 ty			


bus=[...                                                       
     1  1.04000    .0000   0.716   0.27	0		0		1	   ;
     2  1.02500   0   1.63    0.067	0		0			2	   ;
     3  1.02500   0   0.85    -0.109	0		0		2	   ;
     4  1.0260    0 	 0		  0		0		0		3	   ;
     5  0.996   0	 0		  0		1.25    0.50 		3	   ;
     6  1.013     0 	 0		  0		0.9   0.3		3	   ;
     7  1.026     0 	 0		  0		0		0		3	   ;
     8  1.016    0	 0		  0		1    0.35  		3	   ;
     9  1.032     0 	 0		  0		0		0			3	   ;
 ];
  
%qmx and qmn set to +999 -999 default; Vmx Vmn set to 1.5 .5 default

%	fb		tb		r			x		 charging	tap		
 line = [...
     1      4      .00000    .05760    .00000      .00         ;
     2      7      .00000    .0625  .00000      .00          ; 
     3      9      .00000    .0586    .00000      .00          ; 
     4      5      .01    .085    .088*2     .00          ; 
     4      6      .017    .09200    .079*2      .00        ; 
     5      7      .03200    .16100    .153*2      .00         ;
     6      9      .0390    .170    .179*2      .00         ; 
     7      8      .0085    .072    .0745*2      .00         ;
     8      9      .0119     .1008   .1045*2      .00        ]; 
  % tap always at from bus
  % tap = 1 and phase = 0 by default
 
 
%  generator data

% machine data
% num bus b-mva x_1  r_a  x_d  x'_d  x"_d T'_do T"_do
%                         x_q  x'_q  x"_q T'_qo T"_qo
%                         H    d_o   d_1 type s s apf rpf
 mac_con = [ ...
1 1 100 0  0   0.146    0.0608  0       8.96    0 ...
                        0.0969  0.0969  0       0.31    0 ...
                        23.64   0.0125  0       1       0   0;% classical model H = 22.64 damp 0
2 2 100 0  0   0.8958   0.1198  0       6.0 	0 ...
                        0.8645   0.1969   0       0.535   0  ...
                        6.4     0.0068  0       2       0   0;%salient pole generator (type 2), no saturation, damping is so low as to be negligible
3 3 100 0  0   1.3125   0.1813  0       5.89    0 ...
                        1.2578  0.25    0       0.6     0 ...
                        3.01    0.0048   0      1       0   0;% classical model H = 3.10
                ]; 
            
%type machine-# T_R  K_A  T_A      T_B  T_C      V_Rmax  V_Rmin  K_E...
%                T_E  E_1  S_E(E_1) E_2  S_E(E_2) K_F     T_F            
exc_con = [...               
    1     1       0    20 0.2     0    0        1.0     -1.0    1 ...
                 0.314 3.1  0.33     2.3  0.10     0.063     0.35;
    1     2       0    20 0.2     0    0        1.0     -1.0    1 ...
                 0.314 3.1  0.33     2.3  0.10     0.063     0.35;
    1     3       0    20 0.2     0    0        1.0     -1.0    1 ...
                 0.314 3.1  0.33     2.3  0.10     0.063     0.35;
];

