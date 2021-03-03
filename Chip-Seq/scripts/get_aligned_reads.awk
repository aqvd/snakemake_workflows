BEGIN{s=0}
/aligned exactly 1 time/{ s+= $1}
/aligned >1 times/{s += $1} 
END{print s}