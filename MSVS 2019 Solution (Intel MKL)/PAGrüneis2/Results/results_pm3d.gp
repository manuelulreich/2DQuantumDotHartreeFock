set term png size 1280,960
#set xlabel "X" 
#set ylabel "Y"

set size square
#set view 60,36,1
#set view equal xyz
#set view 0,0,1.0
set view map
#center = (0,0)
#set size 1.2,1.2
set parametric
set angles radians
set mapping cylindrical
#set zrange[0:.01]
unset key
set pm3d depthorder interpolate 1,1

set lmargin at screen 0.05;
set rmargin at screen 0.9;
set bmargin at screen 0.1;
set tmargin at screen 0.95;

prefix1 = "Psi"
prefix2 = "PsiSquared"
suffix = ".png"
file = "d_80x160_dr_0.062500_eigenvecs_200.txt"

#test plot
#set out sprintf("%s%i%s",prefix2,2,suffix)
#sp file u 2:(($4)**2):1 with pm3d 

do for [i=3:102] {
#printing psi
	#set zlabel "Psi(x,y)" rotate by 90 offset -2,0,0
	#set out sprintf("%s%i%s",prefix1,i-2,suffix)
	#sp file u 2:i:1 with pm3d
#printing psi squared
	#set zlabel "Psi(x,y)²" rotate by 90 offset -2,0,0
	set out sprintf("%s%i%s",prefix2,i-2,suffix)
	sp file u 2:(column(i)**2):1 with pm3d 
}

unset out