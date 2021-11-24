sets i "rows" /1/;
sets j "columns" /1*3/;
$call csv2gdx Naive.csv id=d index=1 values=2..lastCol useHeader=y
$gdxIn Naive.gdx
*$load i=dim1
*$load j=dim2
parameter d(i,j) ;
$load d
$gdxIn
display i,j,d;

