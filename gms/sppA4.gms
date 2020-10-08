$ontext
    sppA4 pooling problem data.
    Author: Mohammed Alfaki, Wed Nov  9 15:36:08 2011
$offtext

$eolcom #

# Declare sets
    set i    / 1*45 /;
    set s(i) / 1*20 /;
    set t(i) /31*45 /;
    set k    / 1*24 /;

alias (i,j);

# The arc unit cost c_{ij}
table c(i,j)
         21      22      23      24      25      26      27      28      29      30      31      32      33      34      35      36      37      38      39      40      41      42      43      44      45
  1   29.00    0.00   29.00    0.00    0.00    0.00   29.00   29.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -13.00    0.00    0.00    0.00  -17.00    0.00  -16.00  -19.00    0.00
  2    0.00    0.00    0.00    0.00   12.00   12.00    0.00   12.00   12.00   12.00  -35.00  -29.00    0.00    0.00    0.00    0.00    0.00    0.00  -35.00    0.00    0.00    0.00    0.00    0.00    0.00
  3    0.00   44.00    0.00   44.00   44.00   44.00    0.00    0.00   44.00    0.00    0.00    0.00   -1.00    1.00   -6.00    0.00    2.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
  4   16.00    0.00   16.00    0.00    0.00    0.00   16.00    0.00    0.00   16.00    0.00    0.00  -29.00  -27.00    0.00    0.00    0.00    0.00    0.00  -26.00    0.00    0.00    0.00  -32.00    0.00
  5    0.00   29.00   29.00   29.00   29.00   29.00    0.00   29.00    0.00   29.00    0.00    0.00    0.00    0.00    0.00  -11.00    0.00    0.00    0.00  -13.00  -17.00    0.00    0.00  -19.00  -15.00
  6    0.00    0.00    0.00    0.00    0.00   27.00    0.00    0.00   27.00    0.00    0.00    0.00  -18.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -19.00  -20.00    0.00    0.00    0.00
  7   20.00   20.00   20.00   20.00   20.00   20.00   20.00    0.00    0.00    0.00    0.00  -21.00    0.00    0.00    0.00    0.00  -22.00  -24.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
  8   25.00    0.00   25.00   25.00   25.00    0.00   25.00   25.00   25.00    0.00    0.00    0.00  -20.00    0.00  -25.00    0.00    0.00  -19.00  -22.00    0.00    0.00  -22.00    0.00    0.00  -19.00
  9    0.00   17.00   17.00   17.00   17.00   17.00   17.00    0.00   17.00   17.00    0.00    0.00    0.00    0.00  -33.00    0.00  -25.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -27.00
 10    0.00   30.00    0.00    0.00    0.00   30.00    0.00   30.00   30.00   30.00    0.00    0.00  -15.00    0.00    0.00    0.00    0.00    0.00    0.00  -12.00    0.00    0.00  -15.00    0.00    0.00
 11    0.00   10.00   10.00   10.00    0.00   10.00   10.00    0.00   10.00   10.00    0.00  -31.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -32.00  -36.00    0.00  -35.00    0.00  -34.00
 12    0.00    0.00   41.00    0.00    0.00   41.00   41.00   41.00   41.00   41.00   -6.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00   -5.00    0.00    0.00    0.00    0.00
 13   31.00   31.00    0.00   31.00   31.00    0.00    0.00   31.00    0.00    0.00    0.00  -10.00  -14.00    0.00  -19.00    0.00  -11.00    0.00    0.00    0.00  -15.00    0.00    0.00  -17.00    0.00
 14   48.00   48.00    0.00   48.00   48.00   48.00   48.00    0.00   48.00   48.00    0.00    0.00    0.00    0.00   -2.00    0.00    0.00    0.00    0.00    6.00    2.00    0.00    0.00    0.00    4.00
 15    0.00   42.00    0.00   42.00    0.00    0.00   42.00    0.00    0.00   42.00    0.00    1.00    0.00   -1.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
 16    0.00    0.00   36.00   36.00   36.00   36.00   36.00    0.00   36.00    0.00    0.00    0.00    0.00   -7.00    0.00   -4.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
 17   23.00    0.00    0.00    0.00    0.00    0.00    0.00   23.00   23.00    0.00    0.00    0.00  -22.00    0.00    0.00    0.00    0.00    0.00  -24.00    0.00    0.00    0.00    0.00    0.00    0.00
 18   26.00    0.00    0.00    0.00   26.00   26.00   26.00   26.00   26.00   26.00    0.00    0.00    0.00    0.00  -24.00    0.00    0.00    0.00  -21.00    0.00    0.00    0.00    0.00    0.00  -18.00
 19    0.00   37.00   37.00   37.00    0.00   37.00    0.00   37.00    0.00    0.00    0.00    0.00   -8.00    0.00    0.00   -3.00   -5.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00
 20    0.00   42.00   42.00    0.00   42.00    0.00   42.00   42.00    0.00   42.00    0.00    0.00    0.00   -1.00    0.00    0.00    0.00    0.00   -5.00    0.00    0.00    0.00   -3.00    0.00   -2.00
 21    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -47.00    0.00    0.00  -43.00    0.00  -40.00    0.00  -44.00    0.00  -42.00    0.00    0.00  -45.00  -48.00    0.00
 22    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -47.00  -41.00    0.00  -43.00    0.00    0.00  -42.00  -44.00    0.00    0.00    0.00  -47.00    0.00  -48.00    0.00
 23    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -47.00    0.00  -45.00    0.00    0.00    0.00    0.00    0.00    0.00  -42.00    0.00    0.00    0.00    0.00  -44.00
 24    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -47.00  -41.00  -45.00  -43.00    0.00    0.00  -42.00    0.00  -47.00    0.00  -46.00    0.00  -45.00  -48.00    0.00
 25    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -41.00  -45.00  -43.00  -50.00    0.00    0.00  -44.00    0.00    0.00    0.00  -47.00  -45.00    0.00    0.00
 26    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -47.00    0.00    0.00    0.00  -50.00    0.00    0.00    0.00  -47.00  -42.00    0.00    0.00    0.00    0.00  -44.00
 27    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -43.00  -50.00  -40.00  -42.00  -44.00    0.00  -42.00  -46.00    0.00  -45.00  -48.00  -44.00
 28    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -45.00  -43.00    0.00  -40.00  -42.00    0.00  -47.00    0.00    0.00    0.00  -45.00    0.00    0.00
 29    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -45.00  -43.00    0.00  -40.00  -42.00    0.00    0.00    0.00    0.00    0.00  -45.00  -48.00    0.00
 30    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00    0.00  -50.00  -40.00  -42.00    0.00  -47.00    0.00    0.00    0.00    0.00    0.00  -44.00 ;

# The adjacency matrix (the arcs set A)
table a(i,j)
     21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45
  1   1   0   1   0   0   0   1   1   0   0   0   0   0   0   0   0   1   0   0   0   1   0   1   1   0
  2   0   0   0   0   1   1   0   1   1   1   1   1   0   0   0   0   0   0   1   0   0   0   0   0   0
  3   0   1   0   1   1   1   0   0   1   0   0   0   1   1   1   0   1   0   0   0   0   0   0   0   0
  4   1   0   1   0   0   0   1   0   0   1   0   0   1   1   0   0   0   0   0   1   0   0   0   1   0
  5   0   1   1   1   1   1   0   1   0   1   0   0   0   0   0   1   0   0   0   1   1   0   0   1   1
  6   0   0   0   0   0   1   0   0   1   0   0   0   1   0   0   0   0   0   0   0   1   1   0   0   0
  7   1   1   1   1   1   1   1   0   0   0   0   1   0   0   0   0   1   1   0   0   0   0   0   0   0
  8   1   0   1   1   1   0   1   1   1   0   0   0   1   0   1   0   0   1   1   0   0   1   0   0   1
  9   0   1   1   1   1   1   1   0   1   1   0   0   0   0   1   0   1   0   0   0   0   0   0   0   1
 10   0   1   0   0   0   1   0   1   1   1   0   0   1   0   0   0   0   0   0   1   0   0   1   0   0
 11   0   1   1   1   0   1   1   0   1   1   0   1   0   0   0   0   0   0   0   1   1   0   1   0   1
 12   0   0   1   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   1   0   0   0   0
 13   1   1   0   1   1   0   0   1   0   0   0   1   1   0   1   0   1   0   0   0   1   0   0   1   0
 14   1   1   0   1   1   1   1   0   1   1   0   0   0   0   1   0   0   0   0   1   1   0   0   0   1
 15   0   1   0   1   0   0   1   0   0   1   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0
 16   0   0   1   1   1   1   1   0   1   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0
 17   1   0   0   0   0   0   0   1   1   0   0   0   1   0   0   0   0   0   1   0   0   0   0   0   0
 18   1   0   0   0   1   1   1   1   1   1   0   0   0   0   1   0   0   0   1   0   0   0   0   0   1
 19   0   1   1   1   0   1   0   1   0   0   0   0   1   0   0   1   1   0   0   0   0   0   0   0   0
 20   0   1   1   0   1   0   1   1   0   1   0   0   0   1   0   0   0   0   1   0   0   0   1   0   1
 21   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   1   0   1   0   1   0   0   1   1   0
 22   0   0   0   0   0   0   0   0   0   0   1   1   0   1   0   0   1   1   0   0   0   1   0   1   0
 23   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   1   0   0   0   0   1
 24   0   0   0   0   0   0   0   0   0   0   1   1   1   1   0   0   1   0   1   0   1   0   1   1   0
 25   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   0   0   1   0   0   0   1   1   0   0
 26   0   0   0   0   0   0   0   0   0   0   1   0   0   0   1   0   0   0   1   1   0   0   0   0   1
 27   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   0   1   1   0   1   1   1
 28   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   1   1   0   1   0   0   0   1   0   0
 29   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   1   1   0   0   0   0   0   1   1   0
 30   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   0   1   0   0   0   0   0   1 ;

# Source qualities/terminal quality upper bounds
table q(i,k)
          1       2       3       4       5       6       7       8       9      10      11      12      13      14      15      16      17      18      19      20      21      22      23      24
  1   50.28   56.45   73.79   35.11   69.04   62.31    4.08   56.95   35.82    7.07   16.44   52.49  -50.28  -56.45  -73.79  -35.11  -69.04  -62.31   -4.08  -56.95  -35.82   -7.07  -16.44  -52.49
  2   36.78   56.89   69.33   35.89   45.75    4.08    0.39   39.59   15.11   48.67   45.60   35.95  -36.78  -56.89  -69.33  -35.89  -45.75   -4.08   -0.39  -39.59  -15.11  -48.67  -45.60  -35.95
  3   15.06   66.64   29.06   45.29   33.56   21.01   69.81   65.61    6.51   71.07   29.31    7.27  -15.06  -66.64  -29.06  -45.29  -33.56  -21.01  -69.81  -65.61   -6.51  -71.07  -29.31   -7.27
  4   58.59   14.15   14.18    3.73   34.68   66.75   21.72   10.32   73.99   12.93   52.60   18.06  -58.59  -14.15  -14.18   -3.73  -34.68  -66.75  -21.72  -10.32  -73.99  -12.93  -52.60  -18.06
  5   52.19   34.01   16.89   55.54   26.93   61.98    2.02   32.45   70.74   68.99   26.31   32.01  -52.19  -34.01  -16.89  -55.54  -26.93  -61.98   -2.02  -32.45  -70.74  -68.99  -26.31  -32.01
  6    9.89   76.51   37.79   53.64   62.99   21.79    0.96    1.41   49.11   58.63   63.36   78.80   -9.89  -76.51  -37.79  -53.64  -62.99  -21.79   -0.96   -1.41  -49.11  -58.63  -63.36  -78.80
  7   38.57   77.75   27.82   48.38   23.20   14.97   28.05   18.51   26.36   18.53    4.76   53.70  -38.57  -77.75  -27.82  -48.38  -23.20  -14.97  -28.05  -18.51  -26.36  -18.53   -4.76  -53.70
  8   20.74   78.31   32.44   25.47   26.80   47.96   75.99   14.26   23.73   46.48   19.91   57.33  -20.74  -78.31  -32.44  -25.47  -26.80  -47.96  -75.99  -14.26  -23.73  -46.48  -19.91  -57.33
  9   42.12   39.34   25.40    3.57   54.19    4.07   63.40   22.87   49.41   54.57   72.29   17.97  -42.12  -39.34  -25.40   -3.57  -54.19   -4.07  -63.40  -22.87  -49.41  -54.57  -72.29  -17.97
 10    3.87   47.75   40.17   64.00   49.79   41.69   53.54   28.74   42.70   33.41    3.60   11.46   -3.87  -47.75  -40.17  -64.00  -49.79  -41.69  -53.54  -28.74  -42.70  -33.41   -3.60  -11.46
 11   69.55   16.56    8.92   36.84    9.33   23.57   48.96   60.73    2.44   15.26   70.71    6.60  -69.55  -16.56   -8.92  -36.84   -9.33  -23.57  -48.96  -60.73   -2.44  -15.26  -70.71   -6.60
 12    7.51   16.06   51.90   45.00   29.40   68.27   16.15   76.79    5.03   59.89   51.15   15.50   -7.51  -16.06  -51.90  -45.00  -29.40  -68.27  -16.15  -76.79   -5.03  -59.89  -51.15  -15.50
 13   31.08   73.01   71.02   11.59   37.22   52.00   53.53   14.98   39.98   50.96   74.92   21.36  -31.08  -73.01  -71.02  -11.59  -37.22  -52.00  -53.53  -14.98  -39.98  -50.96  -74.92  -21.36
 14   15.09   38.56   57.54   28.53    3.46   66.32   13.61   52.20    8.27   77.19    2.98   36.57  -15.09  -38.56  -57.54  -28.53   -3.46  -66.32  -13.61  -52.20   -8.27  -77.19   -2.98  -36.57
 15   25.12   57.25   70.90   55.46   56.12   51.22   37.44   77.88   61.65   27.06   32.58    2.81  -25.12  -57.25  -70.90  -55.46  -56.12  -51.22  -37.44  -77.88  -61.65  -27.06  -32.58   -2.81
 16   15.67   25.20   26.68   56.52   53.52   61.43    6.85   69.17   42.71   70.40   24.25   54.53  -15.67  -25.20  -26.68  -56.52  -53.52  -61.43   -6.85  -69.17  -42.71  -70.40  -24.25  -54.53
 17   62.17   39.27   50.84   65.75   63.16   12.18   56.85   63.02   53.08   48.32   49.66   66.85  -62.17  -39.27  -50.84  -65.75  -63.16  -12.18  -56.85  -63.02  -53.08  -48.32  -49.66  -66.85
 18   24.36   15.51   75.46   37.35   22.97    3.40   44.55   38.90   61.59   13.30   68.11   58.22  -24.36  -15.51  -75.46  -37.35  -22.97   -3.40  -44.55  -38.90  -61.59  -13.30  -68.11  -58.22
 19   21.45   73.54   28.75   27.03   12.77   11.54   48.87   58.07   21.10   15.89   64.56    4.02  -21.45  -73.54  -28.75  -27.03  -12.77  -11.54  -48.87  -58.07  -21.10  -15.89  -64.56   -4.02
 20    4.92   38.27    8.58   13.42   26.95   28.06   22.31   78.14   38.65   18.40   56.32   66.77   -4.92  -38.27   -8.58  -13.42  -26.95  -28.06  -22.31  -78.14  -38.65  -18.40  -56.32  -66.77
 31   21.12   31.86   34.47   67.59   51.90   66.73   36.46   83.07   90.59   91.77   56.33   25.87   -1.34   -6.70   -5.41   -7.84  -13.14  -18.46   -4.18  -11.83   -8.47   -4.93   -2.01   -8.91
 32   64.32   94.00   23.43   64.42   49.96   98.07   52.00   92.60   69.49   45.49   41.50   55.34  -14.25   -7.14  -19.17   -7.04   -1.97   -4.60  -15.45   -8.89   -6.86  -10.25   -3.97  -11.97
 33   72.22   91.98   96.26   36.34   87.39   71.94   89.95   90.92   22.53   32.66   59.47   85.27   -1.58  -16.14   -3.69   -1.44  -15.03  -12.69   -7.41   -4.68  -14.61  -11.92  -14.62  -14.62
 34   20.08   50.80   74.45   61.36   39.39   95.94   64.83   98.83   43.42   85.67   31.21   94.68   -4.65   -8.73   -9.41  -14.37  -11.07  -10.27  -17.94   -5.25  -14.28  -19.53   -7.93  -13.87
 35   43.21   32.53   75.13   65.21   81.61   42.86   93.05   58.28   57.96   62.65   76.21   86.96   -1.60  -18.78  -12.00  -18.50   -9.37   -3.78   -1.49  -10.68   -1.72  -19.86  -18.04  -13.28
 36   91.20   65.41   43.18   26.66   81.75   36.44   54.96   97.43   38.09   38.77   30.11   34.46  -19.80   -4.33   -4.09  -11.59   -8.25   -0.77   -3.77   -7.56  -14.27  -13.50  -13.36   -1.96
 37   45.29   64.60   73.45   34.38   39.90   38.46   36.57   28.56   32.38   41.34   96.89   82.59   -3.89  -14.91   -1.74   -2.12  -17.49  -14.57  -19.81   -8.08   -4.90  -11.61   -3.66  -16.08
 38   47.78   29.11   87.65   49.80   84.72   81.76   41.18   45.35   96.11   65.42   28.28   74.49   -5.70   -7.79  -15.34  -18.38  -13.13   -7.54   -1.60  -12.81   -4.79  -15.13  -12.69   -8.05
 39   53.39   78.40   47.49   35.72   74.86   49.46   37.19   67.62   52.62   26.71   64.16   51.52   -8.64  -14.45   -9.15  -11.38   -1.07   -0.46   -3.71  -19.99   -1.92   -5.21  -11.26   -4.96
 40   98.69   85.16   74.55   44.41   76.44   36.78   22.71   93.54   58.60   79.49   22.46   23.48  -19.91  -18.78   -4.29   -4.98  -15.94  -12.80   -0.18   -5.35   -9.27   -1.70  -16.28  -13.96
 41   57.15   38.84   33.09   25.04   85.76   71.30   52.46   44.66   34.51   22.58   94.19   70.89   -7.14  -15.65  -16.84   -0.49   -6.59  -14.94   -2.86   -6.16   -2.10   -6.44  -11.72  -10.67
 42   34.18   22.58   21.34   36.23   87.72   65.23   99.41   78.31   97.77   43.14   96.39   90.88   -5.87   -7.96   -2.50   -8.13   -1.73  -11.56   -7.54  -18.26  -18.99  -16.11  -19.13  -11.05
 43   57.61   41.23   52.85   86.47   79.05   38.73   23.44   46.66   94.20   86.24   40.44   60.79   -8.55   -3.32   -7.18   -5.53   -0.30  -11.05   -1.13   -6.02   -4.18   -8.52   -4.63   -6.06
 44   20.83   53.73   35.88   48.62   65.36   50.71   21.08   38.23   59.11   83.54   39.94   28.87  -16.76   -6.51   -0.29  -10.55   -4.05   -2.96   -4.07  -15.64  -17.37  -13.65   -3.26  -13.48
 45   55.09   43.52   36.63   43.09   54.14   57.98   41.59   26.97   45.38   93.70   65.36   48.75   -4.15  -14.15  -12.96   -1.15   -4.81   -2.31  -12.78  -19.71  -10.73  -12.13   -8.57   -7.42 ;

# Node capacity lower bound
parameter bl(i) /  1     0.00
                   2     0.00
                   3     0.00
                   4     0.00
                   5     0.00
                   6     0.00
                   7     0.00
                   8     0.00
                   9     0.00
                  10     0.00
                  11     0.00
                  12     0.00
                  13     0.00
                  14     0.00
                  15     0.00
                  16     0.00
                  17     0.00
                  18     0.00
                  19     0.00
                  20     0.00
                  31     0.00
                  32     0.00
                  33     0.00
                  34     0.00
                  35     0.00
                  36     0.00
                  37     0.00
                  38     0.00
                  39     0.00
                  40     0.00
                  41     0.00
                  42     0.00
                  43     0.00
                  44     0.00
                  45     0.00 / ;

# Node capacity upper bound
parameter bu(i) /  1   222.00
                   2    32.00
                   3   281.00
                   4   205.00
                   5    94.00
                   6   150.00
                   7    78.00
                   8   234.00
                   9   233.00
                  10   280.00
                  11   288.00
                  12    93.00
                  13    23.00
                  14    49.00
                  15   226.00
                  16    12.00
                  17   237.00
                  18    34.00
                  19   281.00
                  20   199.00
                  21    90.00
                  22    59.00
                  23    78.00
                  24   123.00
                  25   111.00
                  26   103.00
                  27   134.00
                  28    53.00
                  29   174.00
                  30   133.00
                  31   128.00
                  32   273.00
                  33   203.00
                  34    46.00
                  35    40.00
                  36   120.00
                  37   211.00
                  38   185.00
                  39    65.00
                  40   222.00
                  41   278.00
                  42   102.00
                  43   157.00
                  44   168.00
                  45   275.00 / ;

$include xmodel.gms