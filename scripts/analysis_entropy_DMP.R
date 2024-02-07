#!/usr/bin/env R

# https://www.yacinemahdid.com/shannon-entropy-from-theory-to-python/
# https://www.yacinemahdid.com/content/images/2022/01/image-8.png


# of gewoon windowed limma probabiliy

a = c(-1,-1,-1,-1,1,1,1,1,1,1,1)
ra = sort(sample(a) - sample(a))



a = c(-1,-1,-1,-1,1,1,1,1,1,1,1)
da = a[1:length(a)-1] - a[2:length(a)]
pa = round(unlist(lapply(da, function(x) {return ( sum(x < ra)/length(ra) )})),4)
pa
-sum(pa * log(pa))



