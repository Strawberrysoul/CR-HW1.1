# Gnuplot script file for plotting data in file "force.dat"
# This file is called   force.p
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set datafile separator ','
set xrange[0:255]
set log y
set yrange[0.01:1000]
set ytics('10^{-2}' 0.01, '10^{-1}' 0.1,'10^{0}' 1, '10^{1}' 10, '10^{2}' 100, '10^{3}' 1000)
plot "C:/Users/Janina/GitHub/CR1/responseCurve.csv" using 1:2 title 'Red' with line lt rgb "red", \
    "C:/Users/Janina/GitHub/CR1/responseCurve.csv" using 1:3 title 'Green' with line lt rgb "green", \
    "C:/Users/Janina/GitHub/CR1/responseCurve.csv" using 1:4 title 'Blue' with line lt rgb "blue"