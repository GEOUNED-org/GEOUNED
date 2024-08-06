Input Test
C   ______ _______  _____      _     _ __   _ _______ ______  
C  |  ____ |______ |     | ___ |     | | \  | |______ |     \ 
C  |_____| |______ |_____|     |_____| |  \_| |______ |_____/
C Version : 1.0.1     22/10/2023
C FreeCAD Version : 0.20.2 
C
C *************************************************************
C Original Step file : inputSTEP/Torus/codo.stp
C
C Creation Date : 2023-10-31 13:37:07.609087
C Solid Cells   : 1
C Total Cells   : 4
C Surfaces      : 11
C Materials     : 0
C
C **************************************************************
1     0      2 4 -3 1
           imp:n=1.0   imp:p=1.0   
C 
C ##########################################################
C              VOID CELLS
C ##########################################################
C 
2     0      5 7 9 -10 -8 -6 (-1:3:-4:-2)
           imp:n=1.0   imp:p=1.0   
           $Automatic Generated Void Cell. Enclosure(-1.0000000000001454, 2.407109860380112, -1.0, 2.407109860380112, -1.3, 1.3)
           $Enclosed cells : (1)
3     0      -11 (-5:6:-7:8:-9:10)
           imp:n=1.0   imp:p=1.0   
           $Graveyard_in
4     0      11
           imp:n=0     imp:p=0     
           $Graveyard
 
C ##########################################################
C                  SURFACE DEFINITION
C ##########################################################
1      PX   1.1102230e-16
2      PY   0.0000000e+00
3      TZ   0.0000000e+00  0.0000000e+00  0.0000000e+00
            1.0000000e+00  3.0000000e-01  3.0000000e-01
4      TZ   0.0000000e+00  0.0000000e+00  0.0000000e+00
            1.0000000e+00  2.0000000e-01  2.0000000e-01
5      PX  -1.0000000e+00
6      PX   2.4071099e+00
7      PY  -1.0000000e+00
8      PY   2.4071099e+00
9      PZ  -1.3000000e+00
10     PZ   1.3000000e+00
11     S   7.0355493e-01  7.0355493e-01  0.0000000e+00  2.7923045e+00
 
C 
MODE P
VOID 
NPS 1e6
PRDMP 2J -1
C SDEF PAR=P X=D1 Y=D2 Z=D3 
C SI1 -1.0000000e+00 2.4071099e+00 
C SI2 -1.0000000e+00 2.4071099e+00 
C SI3 -1.3000000e+00 1.3000000e+00 
C SP1 0  1 
C SP2 0  1 
C SP3 0  1 
SDEF PAR=P NRM=-1 SUR=11 WGT=2.4494886e+01 DIR=d1
SI1 0 1
SP1 -21 1
F4:P  1 
SD4   2.4674011e-01 