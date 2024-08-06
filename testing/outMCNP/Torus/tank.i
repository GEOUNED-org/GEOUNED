Input Test
C   ______ _______  _____      _     _ __   _ _______ ______  
C  |  ____ |______ |     | ___ |     | | \  | |______ |     \ 
C  |_____| |______ |_____|     |_____| |  \_| |______ |_____/
C Version : 1.0.1     22/10/2023
C FreeCAD Version : 0.20.2 
C
C *************************************************************
C Original Step file : inputSTEP/Torus/tank.stp
C
C Creation Date : 2023-10-31 13:37:01.745063
C Solid Cells   : 12
C Total Cells   : 15
C Surfaces      : 34
C Materials     : 0
C
C **************************************************************
1     0      -18 4 -3 -2
           imp:n=1.0   imp:p=1.0   
2     0      -7 19 5 6
           imp:n=1.0   imp:p=1.0   
3     0      2 4 20 -9 -8 1:-1 2 4 21 -9 -8
           imp:n=1.0   imp:p=1.0   
4     0      -8 -4
           imp:n=1.0   imp:p=1.0   
5     0      6 10 -11 22 -5 1:-1 6 10 -11 23 -5
           imp:n=1.0   imp:p=1.0   
6     0      -6 10
           imp:n=1.0   imp:p=1.0   
7     0      -14 -12 13 15
           imp:n=1.0   imp:p=1.0   
8     0      -15 -12 13
           imp:n=1.0   imp:p=1.0   
9     0      8 16 -13 -9 24 1:-1 8 16 -13 -9 24
           imp:n=1.0   imp:p=1.0   
10    0      12 17 -11 -10 25 1:-1 12 17 -11 -10 25
           imp:n=1.0   imp:p=1.0   
11    0      -16 8 -13 1:-26 8 16 -13 1:-16 -1 -13 8:-1 -26 8 16 -13
           imp:n=1.0   imp:p=1.0   
12    0      12 -10 1 -17:-27 12 17 -10 1:-1 12 -10 -17:-1 -27 12 17 -10
           imp:n=1.0   imp:p=1.0   
C 
C ##########################################################
C              VOID CELLS
C ##########################################################
C 
13    0      32 28 30 -31 -29 -33 (18:3:-4:2) (-19:-6:7:-5) (-4:9:8:-2:(-1:-20) 
           (1:-21)) (8:4) (-6:11:-10:5:(-1:-22) (1:-23)) (-10:6) (12:14:-15:-13)
            (12:15:-13) (-8:13:-24:-16:9:1 -1) (10:-12:-25:-17:11:1 -1) (13:-8:
           (16:-1) (-16:26:-1) (16:1) (-16:26:1)) (-12:10:(-1:17) (-1:-17:27) 
           (1:17) (1:-17:27))
           imp:n=1.0   imp:p=1.0   
           $Automatic Generated Void Cell. Enclosure(897.936106449945, 975.703560470437, -1237.0176300593162, -1159.250176038824, -10.500000000000098, 163.49999999999977)
           $Enclosed cells : (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
14    0      -34 (-28:29:-30:31:-32:33)
           imp:n=1.0   imp:p=1.0   
           $Graveyard_in
15    0      34
           imp:n=0     imp:p=0     
           $Graveyard
 
C ##########################################################
C                  SURFACE DEFINITION
C ##########################################################
1      PY  -1.1981339e+03
2      PZ  -3.0660944e+00
3      S   9.3681983e+02 -1.1981339e+03  4.5080000e+01  5.4580000e+01
4      S   9.3681983e+02 -1.1981339e+03  4.5080000e+01  5.2080000e+01
5      PZ   1.5606609e+02
6      S   9.3681983e+02 -1.1981339e+03  1.0792000e+02  5.2080000e+01
7      S   9.3681983e+02 -1.1981339e+03  1.0792000e+02  5.4580000e+01
8      PZ  -8.6079514e-01
9      TZ   9.3681983e+02 -1.1981339e+03  1.2428338e+01
            1.7435000e+01  1.7565000e+01  1.7565000e+01
10     PZ   1.5386080e+02
11     TZ   9.3681983e+02 -1.1981339e+03  1.4057166e+02
            1.7435000e+01  1.7565000e+01  1.7565000e+01
12     PZ   1.4057166e+02
13     PZ   1.2428338e+01
14     C/Z    936.819833 -1198.133903    35.000000
15     C/Z    936.819833 -1198.133903    32.500000
16     TZ   9.3681983e+02 -1.1981339e+03  1.2428338e+01
            1.7435000e+01  1.5065000e+01  1.5065000e+01
17     TZ   9.3681983e+02 -1.1981339e+03  1.4057166e+02
            1.7435000e+01  1.5065000e+01  1.5065000e+01
18     PZ   4.5080000e+01
19     PZ   1.0792000e+02
20     P    8.5141733e-14  2.8883683e-01 -9.5737834e-01 -3.8922381e+02
21     P    8.1439919e-14 -2.8883683e-01 -9.5737834e-01  3.0290658e+02
22     P    0.0000000e+00  2.8883674e-01  9.5737837e-01 -2.4274482e+02
23     P    0.0000000e+00 -2.8883674e-01  9.5737837e-01  4.4938536e+02
24     K/Z   9.368198e+02 -1.198134e+03 -4.176875e+01     0.359596 1
25     K/Z   9.368198e+02 -1.198134e+03  1.947687e+02     0.359597 -1
26     K/Z   9.368198e+02 -1.198134e+03  1.638056e+01     0.359596 -1
27     K/Z   9.368198e+02 -1.198134e+03  1.366194e+02     0.359597 1
28     PX   8.9793611e+02
29     PX   9.7570356e+02
30     PY  -1.2370176e+03
31     PY  -1.1592502e+03
32     PZ  -1.0500000e+01
33     PZ   1.6350000e+02
34     S   9.3681983e+02 -1.1981339e+03  7.6500000e+01  1.0498019e+02
 
C 
MODE P
VOID 
NPS 1e6
PRDMP 2J -1
C SDEF PAR=P X=D1 Y=D2 Z=D3 
C SI1 8.9793611e+02 9.7570356e+02 
C SI2 -1.2370176e+03 -1.1592502e+03 
C SI3 -1.0500000e+01 1.6350000e+02 
C SP1 0  1 
C SP2 0  1 
C SP3 0  1 
SDEF PAR=P NRM=-1 SUR=34 WGT=3.4622994e+04 DIR=d1
SI1 0 1
SP1 -21 1
F4:P  1 2 3 4 5 6 7 8 9 10 11 12 
SD4   4.3507684e+03 4.3507684e+03 1.7674185e+03 5.9242839e+03 1.7674142e+03 
      5.9242727e+03 6.7934383e+04 4.2521892e+05 7.6716880e+03 7.6716924e+03 
      3.8303335e+04 3.8303347e+04 