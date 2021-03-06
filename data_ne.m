% data for New England Test case

disp('Input Case: New England data')
disp('------------------------------------');
global bus line mac_con pss_con exc_con
% -----------
% bus data
%(1.bus#)(2. volt V )( 3.ang in deg  )( 4. pgen )( 5. qgen )
%(6. pload )(7. qload )(8. bus type)
bus=[
   1   1.048  -9.43   0.0000   0.0000   0.0000   0.0000  3 ;
   2   1.0505 -6.89   0.0000   0.0000   0.0000   0.0000  3 ;
   3   1.0341 -9.73   0.0000   0.0000   3.2200   0.0240  3 ;
   4   1.0116 -10.53  0.0000   0.0000   5.0000   1.8400  3 ;
   5   1.0165 -9.38   0.0000   0.0000   0.0000   0.0000  3 ;
   6   1.0172 -8.68   0.0000   0.0000   0.0000   0.0000  3 ;
   7   1.0067 -10.84  0.0000   0.0000   2.3380   8.4000  3 ;
   8   1.0057 -11.34  0.0000   0.0000   5.2200   1.7600  3 ;
   9   1.0322 -11.15  0.0000   0.0000   0.0000   0.0000  3 ;
  10   1.0235 -6.31   0.0000   0.0000   0.0000   0.0000  3 ;
  11   1.0201 -7.12   0.0000   0.0000   0.0000   0.0000  3 ;
  12   1.0072 -7.14   0.0000   0.0000   0.0850   0.8800  3 ;
  13   1.0207 -7.02   0.0000   0.0000   0.0000   0.0000  3 ;
  14   1.0181 -8.66   0.0000   0.0000   0.0000   0.0000  3 ;
  15   1.0194 -9.06   0.0000   0.0000   3.2000   1.5300  3 ;
  16   1.0346 -7.66   0.0000   0.0000   3.2940   3.2300  3 ;
  17   1.0365 -8.65   0.0000   0.0000   0.0000   0.0000  3 ;
  18   1.0343 -9.49   0.0000   0.0000   1.5800   0.3000  3 ;
  19   1.0509 -3.04   0.0000   0.0000   0.0000   0.0000  3 ;
  20   0.9914 -4.45   0.0000   0.0000   6.8000   1.0300  3 ;
  21   1.0337 -5.26   0.0000   0.0000   2.7400   1.1500  3 ;
  22   1.0509 -0.82   0.0000   0.0000   0.0000   0.0000  3 ;
  23   1.0459 -1.02   0.0000   0.0000   2.4750   0.8460  3 ;
  24   1.0399 -7.54   0.0000   0.0000   3.0860   -0.922  3 ;
  25   1.0587 -5.51   0.0000   0.0000   2.2400   0.4720  3 ;
  26   1.0536 -6.77   0.0000   0.0000   1.3900   0.1700  3 ;
  27   1.0399 -8.78   0.0000   0.0000   2.8100   0.7550  3 ;
  28   1.0509 -3.27   0.0000   0.0000   2.0600   0.2760  3 ;
  29   1.0505 -0.51   0.0000   0.0000   2.8350   1.2690  3 ;
  30   1.0475 -4.47   2.5000   1.3621   0.0000   0.0000  2 ;
  31   1.0400  0.00   5.7293   1.7036   0.0920   0.0460  2 ;
  32   0.9831  1.63   6.5000   1.7590   0.0000   0.0000  2 ;
  33   0.9972  2.18   6.3200   1.0335   0.0000   0.0000  2 ;
  34   1.0123  0.74   5.0800   1.6400   0.0000   0.0000  2 ;
  35   1.0493  4.14   6.5000   2.0884   0.0000   0.0000  2 ;
  36   1.0635  6.83   5.6000   0.9688   0.0000   0.0000  2 ;
  37   1.0278  1.26   5.4000  -0.0444   0.0000   0.0000  2 ;
  38   1.0265  6.55   8.3000   0.1939   0.0000   0.0000  2 ;
  39   1.0300 -10.96  10.000   0.6846   11.040   2.5000  1 ];

% -----------
% line data 
%1.b#1 2.b#2( 3.res  R )( 4. rea X  )(5.charg  B)(6.tap ratio )
line = [
   1   2 0.00350  0.04110 0.69870 0.0000 ;
   1  39 0.00100  0.02500 0.75000 0.0000 ;
   2   3 0.00130  0.01510 0.25720 0.0000 ;
   2  25 0.00700  0.00860 0.14600 0.0000 ;
   2  30 0.00000  0.01810 0.00000 1.0250 ;
   3   4 0.00130  0.02130 0.22140 0.0000 ;
   3  18 0.00110  0.01330 0.21380 0.0000 ;
   4   5 0.00080  0.01280 0.13420 0.0000 ;
   4  14 0.00080  0.01290 0.13820 0.0000 ;
   5   8 0.00080  0.01120 0.14760 0.0000 ;
   6   5 0.00020  0.00260 0.04340 0.0000 ;
   6   7 0.00060  0.00920 0.11300 0.0000 ;
   6  11 0.00070  0.00820 0.13890 0.0000 ;
   7   8 0.00040  0.00460 0.07800 0.0000 ;
   8   9 0.00230  0.03630 0.38040 0.0000 ;
   9  39 0.00100  0.02500 1.20000 0.0000 ;
  10  11 0.00040  0.00430 0.07290 0.0000 ;
  10  13 0.00040  0.00430 0.07290 0.0000 ;
  10  32 0.00000  0.02000 0.00000 1.0700 ;
  12  11 0.00160  0.04350 0.00000 1.0060 ;
  12  13 0.00160  0.04350 0.00000 1.0060 ;
  13  14 0.00090  0.01010 0.17230 0.0000 ;
  14  15 0.00180  0.02170 0.36600 0.0000 ;
  15  16 0.00090  0.00940 0.17100 0.0000 ;
  16  17 0.00070  0.00890 0.13420 0.0000 ;
  16  19 0.00160  0.01950 0.30400 0.0000 ;
  16  21 0.00080  0.01350 0.25480 0.0000 ;
  16  24 0.00030  0.00590 0.06800 0.0000 ;
  17  18 0.00070  0.00820 0.13190 0.0000 ;
  17  27 0.00130  0.01730 0.32160 0.0000 ;
  19  33 0.00070  0.01420 0.00000 1.0700 ;
  19  20 0.00070  0.01380 0.00000 1.0600 ;
  20  34 0.00090  0.01800 0.00000 1.0090 ;
  21  22 0.00080  0.01400 0.25650 0.0000 ;
  22  23 0.00060  0.00960 0.18460 0.0000 ;
  22  35 0.00000  0.01430 0.00000 1.0250 ;
  23  24 0.00220  0.03500 0.36100 0.0000 ;
  23  36 0.00050  0.02720 0.00000 1.0000 ;
  25  26 0.00320  0.03230 0.51300 0.0000 ;
  25  37 0.00060  0.02320 0.00000 1.0250 ;
  26  27 0.00140  0.01470 0.23960 0.0000 ;
  26  28 0.00430  0.04740 0.78020 0.0000 ;
  26  29 0.00570  0.06250 1.02900 0.0000 ;
  28  29 0.00140  0.01510 0.24900 0.0000 ;
  29  38 0.00080  0.01560 0.00000 1.0250 ;
  31   6 0.00000  0.02500 0.00000 1.0000 ];

% -----------
% machine data
% num bus b-mva x_1  r_a  x_d  x'_d  x"_d T'_do T"_do
%                         x_q  x'_q  x"_q T'_qo T"_qo
%                         H    d_o   d_1 types  apf rpf
mac_con = [
  1  30  1000.0 0.130 0  1.000 0.31  0.31 10.20 0.00 ...
                         0.690 0.69  0.69  1.500   0 ...
                         4.200 0     0.00 2 0 0;
  2  31  1000.0 0.350 0  2.950 0.7   0.7 6.560 0.00 ...
                         2.820 1.7   1.7 1.500 0.00 ...
                         3.030 0     0.00 3 0 0;
  3  32  1000.0 0.30  0  2.500 0.530 0.530 5.700 0.0 ... 
                         2.370 0.88  0.88   1.500 0.0 ...
                         3.580 0      0.00 3 0 0;
  4  33  1000.0 0.295 0  2.620 0.440 0.440 5.690 0.0 ...
                         2.580 1.66  1.66  1.500 0.0 ...
                         2.860 0     0.00 3 0 0;
  5  34  1000.0 0.540 0  6.700 1.320 1.32 5.400 0.0 ...
                         6.200 1.660 1.66 0.440 0.0 ...
                         2.600 0     0.00 3 0 0;
  6  35  1000.0 0.224 0  2.540 0.500  0.5  7.300 0.0 ...
                         2.410 0.810 0.81  0.400 0.0 ... 
                         3.480 0     0.00 3 0 0;
  7  36  1000.0 0.322 0  2.900 0.490 0.49 5.660 0.0 ...
                         2.800 1.86  1.86 1.500 0.0 ...
                         2.640 0     0.00 3 0 0;
  8  37  1000.0 0.280 0  2.900 0.570 0.570 6.700 0.0 ...
                         2.800 0.91 0.91 0.410 0.0 ...
                         2.430 0     0.00 3 0 0;
  9  38  1000.0 0.298 0  2.106 0.570 0.570 4.790 0.0 ...
                         2.050 0.590 0.590 1.960 0.0 ...
                         3.450 0     0.00 3 0 0;
  10 39  1000.0 0.030 0  0.200 0.060 0.06  7.000 0 ...
                         0.190 0.080 0.08  0.700 0 ...
                         500 0 0.00 3 0 0;
];

% -----------
% exciter data data from IEEE recommended practice for excitation system
% model for power system stability study
% DC1A exciter 
% type machine-# T_R  K_A  T_A      T_B  T_C      V_Rmax  V_Rmin  K_E...
%                T_E  E_1  S_E(E_1) E_2  S_E(E_2) K_F     T_F
exc_con = [
   1     1       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;
   1     2       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;
   1     3       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;
   1     4       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;
   1     5       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;
   1     6       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;
   1     7       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;    
   1     8       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;   
   1     9       0    46.0 0.06     0    0        1.0     -1.0    0 ...
                 0.46 3.1  0.33     2.3  0.10     0.1     1.0;                 
  ];




% -----------
% pss data type PSS1A with speed input
% column 1 type 1 for speed input 
% column 2 generator number
% column 3 Gpss pss gain
% column 4 Tw
% column 5 Tn1
% column 6 Td1
% column 7 Tn2
% column 8 Td2
% column 9 ymax
% column 10 ymin
pss_con=[
    1  1  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    1  2  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    1  3  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    1  4  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    1  5  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    1  6  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    1  7  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    1  8  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    1  9  3.15  10.0  0.76  0.1  0.76  0.1  0.09  -0.09;
    ];




