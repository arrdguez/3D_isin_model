 
set terminal pdfcairo
set output './plot_total_m.pdf'

#set format y '%.2t*10^%+03T'

set grid
plot './Data.dat'  u 1:($4):5 w errorbar pt 7 ps 0.5, '' u 1:($4) w lp pt 7 ps 0.5


Ti = 1
Tf = 10
a = 0.2
cont=(Tf-Ti)/a
#iniT + i * Tsteps

do for [i=0:cont] {
  c = Ti + i * a
  set title sprintf('temp=%f',c)
  unset xlabel  

  set output sprintf('test%f.pdf',c)
  set pointsize 0.3
  
  set multiplot layout 2,1 rowsfirst

    @TMARGIN; @LMARGIN; @RMARGIN 
    

    set ylabel "Energy/1e18"
    p 'liveData.dat' u (abs($2-(c))<0.01?$1:0/0):3 w p pt 7 title '',\
    '' u (abs($2-(c))<0.01?$1:0/0):4 w lp pt 7 t ""
  
    set xlabel "Time"
    set ylabel "Magnetization"
    unset title
    @BMARGIN; @LMARGIN; @RMARGIN
    p 'liveData.dat' u (abs($2-(c))<0.01?$1:0/0):6 w p pt 7 title '',\
    '' u (abs($2-(c))<0.01?$1:0/0):7 w lp pt 7 t ""

  unset multiplot
  
}
