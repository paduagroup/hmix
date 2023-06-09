#! /usr/local/bin/gnuplot

# set size ratio 1

set xlabel 'x_2'
set xtics 0.2

set ylabel 'h_i / kJ mol^{-1}'
set ytics 2.0
set yrange[-11.0:7.0]

plot 0 w l lt 0 t '',\
     'dioxane_water_h1.out' u 1:2 lt 1 t 'h_1',\
     'dioxane_water_h2.out' u 1:2 lt 2 t 'h_2',\
     'dioxane_water_hrk.mat' u 1:($2/1000) w l lt 3 t 'h_1 Mat',\
     'dioxane_water_hrk.mat' u 1:($3/1000) w l lt 4 t 'h_2 Mat',\
     'dioxane_water_hrk.out' u 1:2 w l lt 1  t 'h_1',\
     'dioxane_water_hrk.out' u 1:3 w l lt 2  t 'h_2',\
     'dioxane_water_hrk.out' u 1:4 w l lt -1 t 'H^E'
     
pause -1
