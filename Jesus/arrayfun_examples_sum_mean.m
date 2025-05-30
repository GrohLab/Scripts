% introduction to anonymous functions

 out = arrayfun(@(var2) arrayfun(@(var1) sum(x(var1,var2,1:2:end),3), 1:size(x,1)), 1:size(x,2),'UniformOutput',false)
 out = arrayfun(@(var2) arrayfun(@(var1) sum(x(var1,var2,1:2:end),3), 1:size(x,1), 'UniformOutput', false), 1:size(x,2), 'UniformOutput', false)
 
 
out = cat(1, out{:});
out

out =

     2     3     4     5     6     7     8     9    10    11
     3     4     5     6     7     8     9    10    11    12
     4     5     6     7     8     9    10    11    12    13
     5     6     7     8     9    10    11    12    13    14
     6     7     8     9    10    11    12    13    14    15
     7     8     9    10    11    12    13    14    15    16
     8     9    10    11    12    13    14    15    16    17
     9    10    11    12    13    14    15    16    17    18
    10    11    12    13    14    15    16    17    18    19
    11    12    13    14    15    16    17    18    19    20
    12    13    14    15    16    17    18    19    20    21
    13    14    15    16    17    18    19    20    21    22

out = arrayfun(@(var2) arrayfun(@(var1) var1 + var2, 1:size(x,1)), 1:size(x,2), 'UniformOutput', false);
out = arrayfun(@(var2) arrayfun(@(var1) x(var1) + x(var2), 1:size(x,1)), 1:size(x,2), 'UniformOutput', false);
 arrayfun(@(var2) arrayfun(@(var1) x(var1) + x(var2), 1:size(x,1)), 1:size(x,2), 'UniformOutput', false);
x

x(:,:,1) =

    43    54    78    52    26    92    18    27    65    46    59    55
    10    66    43    95    23     1    73    77    68    67    55    65
    60    41    10    64    67    47    48    19    64    78    87    55
    48    82    27    96    85    43    16    29    95    36    27    73
    70    72    16    25    35    47    35    10    21    67    32    53
    70    97    29    68    79    78    61    58    71    42    12   100
    64    54    45    29    68    33    20    69    24    85    94    22
     4    33    53    68     1    79    74    55    12    84    65    11
     7    11    46    70    61    48    25    43    61    26    48    11
    32    62    88     7    39     4    92    65    46    62    64     7


x(:,:,2) =

    41    70    35    74    83    80    52    54    86    62    74    77
    45    10    15    40    43    95    89     9    57    99    59    59
    37    53    59    69    89    33    59    12    93    53    25    93
    77    54    27    71    40    68    16    14    70    48    67    59
    63    87     5    45    77    44    20    68    59    81     9     2
    78    49    76     2    40    84    41    50    82    23    63    13
    94    40    25    34    81    77    75    19    88    50    67    87
    98    68    45    43    76    17    83    50    99    91    73    49
    20    75    69    28    38    87    79    15     1    58    90    85
    14    53    36    20    22    99    32     6    87    85    99    21


x(:,:,3) =

    56    15    13    95    74    14    56    99    18    92    90    46
    63    19    50     9     7     4    86    54    36    11     8    11
     4     5    86    11    87    94    35    71     6    75    25   100
    62    64    88    15    94    31    45   100    53    74     6    34
    37    29    28    17    99    30     6    29    34    57    45    30
     5    54    21    63    86    34    18    42    18    19     2     7
    49    70    57    58    79    47    67    47    21    60    90    30
    20    50    65     6    52    65    34    77    91    30    20     5
    13    54    42    94    18     3    90    82    68    14    10    51
    21    45    21    73    40    85    12    11    47    22    31    77

arrayfun(@(var2) arrayfun(@(var1) x(var1) + x(var2), 1:size(x,1)), 1:size(x,2), 'UniformOutput', false);
out = arrayfun(@(var2) arrayfun(@(var1) x(var1) + x(var2), 1:size(x,1)), 1:size(x,2), 'UniformOutput', false);
out = arrayfun(@(var2) arrayfun(@(var1) x(var1) + x(var2), 1:size(x,1)), 1:size(x,2), 'UniformOutput', false);
out = arrayfun(@(var2) arrayfun(@(var1) x(var1) + x(var2), 1:size(x,1)), 1:size(x,2), 'UniformOutput', false);
out = arrayfun(@(var2) arrayfun(@(var1) x(var1) + x(var2), 1:size(x,1)), 1:size(x,2), 'UniformOutput', false)

out =

  1�12 cell array

  Columns 1 through 7

    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}

  Columns 8 through 12

    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}

out = arrayfun(@(var2) arrayfun(@(var1) mean(x(var1,var2,:),3), 1:size(x,1)), 1:size(x,2), 'UniformOutput', false)

out =

  1�12 cell array

  Columns 1 through 7

    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}

  Columns 8 through 12

    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}

out = cat(1, out{:});
out = cat(1, out{:})
Brace indexing is not supported for variables of this type.
 
out

out =

   46.6667   39.3333   33.6667   62.3333   56.6667   51.0000   69.0000   40.6667   13.3333   22.3333
   46.3333   31.6667   33.0000   66.6667   62.6667   66.6667   54.6667   50.3333   46.6667   53.3333
   42.0000   36.0000   51.6667   47.3333   16.3333   42.0000   42.3333   54.3333   52.3333   48.3333
   73.6667   48.0000   48.0000   60.6667   29.0000   44.3333   40.3333   39.0000   64.0000   33.3333
   61.0000   24.3333   81.0000   73.0000   70.3333   68.3333   76.0000   43.0000   39.0000   33.6667
   62.0000   33.3333   58.0000   47.3333   40.3333   65.3333   52.3333   53.6667   46.0000   62.6667
   42.0000   82.6667   47.3333   25.6667   20.3333   40.0000   54.0000   63.6667   64.6667   45.3333
   60.0000   46.6667   34.0000   47.6667   35.6667   50.0000   45.0000   60.6667   46.6667   27.3333
   56.3333   53.6667   54.3333   72.6667   38.0000   57.0000   44.3333   67.3333   43.3333   60.0000
   66.6667   59.0000   68.6667   52.6667   68.3333   28.0000   65.0000   68.3333   32.6667   56.3333
   74.3333   40.6667   45.6667   33.3333   28.6667   25.6667   83.6667   52.6667   49.3333   64.6667
   59.3333   45.0000   82.6667   55.3333   28.3333   40.0000   46.3333   21.6667   49.0000   35.0000

x

x(:,:,1) =

    43    54    78    52    26    92    18    27    65    46    59    55
    10    66    43    95    23     1    73    77    68    67    55    65
    60    41    10    64    67    47    48    19    64    78    87    55
    48    82    27    96    85    43    16    29    95    36    27    73
    70    72    16    25    35    47    35    10    21    67    32    53
    70    97    29    68    79    78    61    58    71    42    12   100
    64    54    45    29    68    33    20    69    24    85    94    22
     4    33    53    68     1    79    74    55    12    84    65    11
     7    11    46    70    61    48    25    43    61    26    48    11
    32    62    88     7    39     4    92    65    46    62    64     7


x(:,:,2) =

    41    70    35    74    83    80    52    54    86    62    74    77
    45    10    15    40    43    95    89     9    57    99    59    59
    37    53    59    69    89    33    59    12    93    53    25    93
    77    54    27    71    40    68    16    14    70    48    67    59
    63    87     5    45    77    44    20    68    59    81     9     2
    78    49    76     2    40    84    41    50    82    23    63    13
    94    40    25    34    81    77    75    19    88    50    67    87
    98    68    45    43    76    17    83    50    99    91    73    49
    20    75    69    28    38    87    79    15     1    58    90    85
    14    53    36    20    22    99    32     6    87    85    99    21


x(:,:,3) =

    56    15    13    95    74    14    56    99    18    92    90    46
    63    19    50     9     7     4    86    54    36    11     8    11
     4     5    86    11    87    94    35    71     6    75    25   100
    62    64    88    15    94    31    45   100    53    74     6    34
    37    29    28    17    99    30     6    29    34    57    45    30
     5    54    21    63    86    34    18    42    18    19     2     7
    49    70    57    58    79    47    67    47    21    60    90    30
    20    50    65     6    52    65    34    77    91    30    20     5
    13    54    42    94    18     3    90    82    68    14    10    51
    21    45    21    73    40    85    12    11    47    22    31    77

x>50

  10�12�3 logical array

ans(:,:,1) =

   0   1   1   1   0   1   0   0   1   0   1   1
   0   1   0   1   0   0   1   1   1   1   1   1
   1   0   0   1   1   0   0   0   1   1   1   1
   0   1   0   1   1   0   0   0   1   0   0   1
   1   1   0   0   0   0   0   0   0   1   0   1
   1   1   0   1   1   1   1   1   1   0   0   1
   1   1   0   0   1   0   0   1   0   1   1   0
   0   0   1   1   0   1   1   1   0   1   1   0
   0   0   0   1   1   0   0   0   1   0   0   0
   0   1   1   0   0   0   1   1   0   1   1   0


ans(:,:,2) =

   0   1   0   1   1   1   1   1   1   1   1   1
   0   0   0   0   0   1   1   0   1   1   1   1
   0   1   1   1   1   0   1   0   1   1   0   1
   1   1   0   1   0   1   0   0   1   0   1   1
   1   1   0   0   1   0   0   1   1   1   0   0
   1   0   1   0   0   1   0   0   1   0   1   0
   1   0   0   0   1   1   1   0   1   0   1   1
   1   1   0   0   1   0   1   0   1   1   1   0
   0   1   1   0   0   1   1   0   0   1   1   1
   0   1   0   0   0   1   0   0   1   1   1   0


ans(:,:,3) =

   1   0   0   1   1   0   1   1   0   1   1   0
   1   0   0   0   0   0   1   1   0   0   0   0
   0   0   1   0   1   1   0   1   0   1   0   1
   1   1   1   0   1   0   0   1   1   1   0   0
   0   0   0   0   1   0   0   0   0   1   0   0
   0   1   0   1   1   0   0   0   0   0   0   0
   0   1   1   1   1   0   1   0   0   1   1   0
   0   0   1   0   1   1   0   1   1   0   0   0
   0   1   0   1   0   0   1   1   1   0   0   1
   0   0   0   1   0   1   0   0   0   0   0   1

arrayfun(@(varEm) varEm>50, x(:,:,:))

  10�12�3 logical array

ans(:,:,1) =

   0   1   1   1   0   1   0   0   1   0   1   1
   0   1   0   1   0   0   1   1   1   1   1   1
   1   0   0   1   1   0   0   0   1   1   1   1
   0   1   0   1   1   0   0   0   1   0   0   1
   1   1   0   0   0   0   0   0   0   1   0   1
   1   1   0   1   1   1   1   1   1   0   0   1
   1   1   0   0   1   0   0   1   0   1   1   0
   0   0   1   1   0   1   1   1   0   1   1   0
   0   0   0   1   1   0   0   0   1   0   0   0
   0   1   1   0   0   0   1   1   0   1   1   0


ans(:,:,2) =

   0   1   0   1   1   1   1   1   1   1   1   1
   0   0   0   0   0   1   1   0   1   1   1   1
   0   1   1   1   1   0   1   0   1   1   0   1
   1   1   0   1   0   1   0   0   1   0   1   1
   1   1   0   0   1   0   0   1   1   1   0   0
   1   0   1   0   0   1   0   0   1   0   1   0
   1   0   0   0   1   1   1   0   1   0   1   1
   1   1   0   0   1   0   1   0   1   1   1   0
   0   1   1   0   0   1   1   0   0   1   1   1
   0   1   0   0   0   1   0   0   1   1   1   0


ans(:,:,3) =

   1   0   0   1   1   0   1   1   0   1   1   0
   1   0   0   0   0   0   1   1   0   0   0   0
   0   0   1   0   1   1   0   1   0   1   0   1
   1   1   1   0   1   0   0   1   1   1   0   0
   0   0   0   0   1   0   0   0   0   1   0   0
   0   1   0   1   1   0   0   0   0   0   0   0
   0   1   1   1   1   0   1   0   0   1   1   0
   0   0   1   0   1   1   0   1   1   0   0   0
   0   1   0   1   0   0   1   1   1   0   0   1
   0   0   0   1   0   1   0   0   0   0   0   1

out = arrayfun(@(var2) arrayfun(@(var1) x(var1,var2,:) > 50, 1:size(x,1)), 1:size(x,2), 'UniformOutput', false)
Error using arrayfun
Non-scalar in Uniform output, at index 1, output 1.
Set 'UniformOutput' to false.

Error in (var2)arrayfun(@(var1)x(var1,var2,:)>50,1:size(x,1))
 
out = arrayfun(@(var2) arrayfun(@(var1) x(var1,var2,:) > 50, 1:size(x,1), 'UniformOutput', false), 1:size(x,2), 'UniformOutput', false)

out =

  1�12 cell array

  Columns 1 through 8

    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}

  Columns 9 through 12

    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}

out = cat(1, out{:})

out =

  12�10 cell array

  Columns 1 through 6

    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}

  Columns 7 through 10

    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}
    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}    {1�1�3 logical}

x

x(:,:,1) =

    43    54    78    52    26    92    18    27    65    46    59    55
    10    66    43    95    23     1    73    77    68    67    55    65
    60    41    10    64    67    47    48    19    64    78    87    55
    48    82    27    96    85    43    16    29    95    36    27    73
    70    72    16    25    35    47    35    10    21    67    32    53
    70    97    29    68    79    78    61    58    71    42    12   100
    64    54    45    29    68    33    20    69    24    85    94    22
     4    33    53    68     1    79    74    55    12    84    65    11
     7    11    46    70    61    48    25    43    61    26    48    11
    32    62    88     7    39     4    92    65    46    62    64     7


x(:,:,2) =

    41    70    35    74    83    80    52    54    86    62    74    77
    45    10    15    40    43    95    89     9    57    99    59    59
    37    53    59    69    89    33    59    12    93    53    25    93
    77    54    27    71    40    68    16    14    70    48    67    59
    63    87     5    45    77    44    20    68    59    81     9     2
    78    49    76     2    40    84    41    50    82    23    63    13
    94    40    25    34    81    77    75    19    88    50    67    87
    98    68    45    43    76    17    83    50    99    91    73    49
    20    75    69    28    38    87    79    15     1    58    90    85
    14    53    36    20    22    99    32     6    87    85    99    21


x(:,:,3) =

    56    15    13    95    74    14    56    99    18    92    90    46
    63    19    50     9     7     4    86    54    36    11     8    11
     4     5    86    11    87    94    35    71     6    75    25   100
    62    64    88    15    94    31    45   100    53    74     6    34
    37    29    28    17    99    30     6    29    34    57    45    30
     5    54    21    63    86    34    18    42    18    19     2     7
    49    70    57    58    79    47    67    47    21    60    90    30
    20    50    65     6    52    65    34    77    91    30    20     5
    13    54    42    94    18     3    90    82    68    14    10    51
    21    45    21    73    40    85    12    11    47    22    31    77

x=randi(100,10,12,5)

x(:,:,1) =

    64    75    70    55    69    75    79    49    12    11    68    10
     9     2    56    49    14    24    37    16    79    94    43     1
     9     5    40    90    73    74    21    79    30    19    46    43
    78    67     7    80    12    98     9    11    61    27    61    66
    91    61    79    74    12    87    78    30    97    80     6    73
    54    53    34     6    65     9    21    24    44    49    32    54
    11    73    61     8    33    37    39    54    70    77    78    11
    83    71    75     9    66    37    56    10    76    40    70    64
    34    79    11    80    75    69    23    41    44    28    13    13
    30    29    13    95    59    60    65    11    66     4    14    14


x(:,:,2) =

    10    56    99    16    60    13    76    17    81     3    76    32
    15    19    18    39    34     3    75    67    75    93    23    82
    17    22    26    17    30    30    75    90    13    66     7    79
    20     8    40    76    46    32    11    52    53    94    77    86
    32    92     8    88    43    66    69    71    33    17    68    51
    32    71    69    36    36    96    47    16    55    93    72    64
    22    56    41    69    56    94    22    96    40    80    65    96
    26    32    99    30    75    46    10    55    42    58    42    45
    90    17    41    54    43    25    83    68    19    45    40     7
    71    63    63    84    43    77    18     4    26    26    82    87


x(:,:,3) =

    64    19    23    63    28    99    14    95    42    36    79    23
    36    73    38     3    25     7    22    68    61    98    70    27
   100    38     9    92    46    94    19    99    76    35     1    68
    23    85    65    81    23     2     5    77    59    89    85    48
    66    74    19    75    81    69    11    34    56    46    93    63
    61    58     5    82    99    79    62    67    59    42    78    24
    39    18    73    39     3    54    94    25    52    22     5    18
    15    96    35    62    54    89    36    30     9    13    38    83
     3    27    67    58     9    90    42    69    72    31    71    77
    43    93    39    54    81    63    99    53   100    73    73    94


x(:,:,4) =

    11    32    51    72    84    50    29    72    51    59    41    94
    19    18    44    62    33    70    24    86    49    68    13    40
    10    34   100    35    56    98    72    29    88    37    27     5
    49    22    82    94    98    33    63    74    36    63    26    35
    20    52    49    13    55    84    60    14    45    82    34    74
    90    91    90    74    34    74    67    84    97     2    16    80
    10    63    14    65    62    96     5    14     5     9    35    55
     5    11    40    84    37     4    35    59    98    98    13    69
    56    40    93    40    76    36    46    37    19    66    89    90
    78     6    92    75    42    67    25    81    67    24    10     6


x(:,:,5) =

    31    29    51     3    95    83    41    39    46    28    58    12
     5    55    65    56    55    85    67    46    21    72    33    82
    20    99    31    31    73    38    94    25    90    29    46    33
    73    72    14    94    58    60    82    79    77    90    72    25
    73    84    48    99     3    88    49    89    89    83    89    35
    88    44    37    29    45    94    76    92    29    40    73    38
    59    48    79    81    65    67    42    56    68    50     2    55
     8    57    79    90    53    21    98    60    67    70    68    57
    93    27    67    60    38    66    99    15    13    84    44    40
    81    75    14    89    94     8    87    90    41    61    44    40

out = arrayfun(@(var2) arrayfun(@(var1) sum(x(var1,var2,1:2:end),3), 1:size(x,1), 'UniformOutput', false), 1:size(x,2), 'UniformOutput', false)

out =

  1�12 cell array

    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}    {1�10 cell}

out = arrayfun(@(var2) arrayfun(@(var1) sum(x(var1,var2,1:2:end),3), 1:size(x,1)), 1:size(x,2))
Error using arrayfun
Non-scalar in Uniform output, at index 1, output 1.
Set 'UniformOutput' to false.
 
out = arrayfun(@(var2) arrayfun(@(var1) sum(x(var1,var2,1:2:end),3), 1:size(x,1)), 1:size(x,2),'UniformOutput',false)

out =

  1�12 cell array

    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}

out = cat(2, out{:})

out =

  Columns 1 through 35

   159    50   129   174   230   203   109   106   130   154   123   130   142   224   219   155   139   224   133   197   144   159    80    86   146    76   213   189   145    66   121   108   213   255   248

  Columns 36 through 70

   117   128   161   198   238   192    94   192    93    96   209   101   173   122   234   257   116   206   160   244   182   158   147   225   131   134   126   134    96   138   159   175   190   164   251

  Columns 71 through 105

   183   130   203   167   153   183   135   100   125   154   100   161   196   197   242   132   190   152   129   207    75   264    83   206   209   131   149   123   143   138   205   146    93   218   188

  Columns 106 through 120

   183    85   176   128   131    45   110   144   139   171   116    84   204   130   148

out = cat(1, out{:})
Brace indexing is not supported for variables of this type.
 
out = arrayfun(@(var2) arrayfun(@(var1) sum(x(var1,var2,1:2:end),3), 1:size(x,1)), 1:size(x,2),'UniformOutput',false)

out =

  1�12 cell array

    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}    {1�10 double}

out = cat(1, out{:})

out =

   159    50   129   174   230   203   109   106   130   154
   123   130   142   224   219   155   139   224   133   197
   144   159    80    86   146    76   213   189   145    66
   121   108   213   255   248   117   128   161   198   238
   192    94   192    93    96   209   101   173   122   234
   257   116   206   160   244   182   158   147   225   131
   134   126   134    96   138   159   175   190   164   251
   183   130   203   167   153   183   135   100   125   154
   100   161   196   197   242   132   190   152   129   207
    75   264    83   206   209   131   149   123   143   138
   205   146    93   218   188   183    85   176   128   131
    45   110   144   139   171   116    84   204   130   148

size(out)

ans =

    12    10

size(x)

ans =

    10    12     5

size(out)

ans =

    12    10

out = arrayfun(@(var2) arrayfun(@(var1) sum(x(var1,var2,1:2:end),3), 1:size(x,1), 'UniformOutput', false), 1:size(x,2), 'UniformOutput', false)