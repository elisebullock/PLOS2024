library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}

f=function(x){
  mix=c(0.8, 0.02, 0.02,0.16)
  timescale=c(4, 30, 60, 1000)
  y=0
  for(i in 1:length(mix)){
  y=y+mix[i]*exp(-x/timescale[i])  
  }
  return(y)
}
g=function(x){0.9*exp(-x/5) + 0.1*exp(-x/(20+x/5))}
tp=seq(0,500,0.1)

d=data.frame(x=tp, y=g(tp))

p=ggplot(d, aes(x, y=y)) + geom_line() +scale_y_log10(minor_breaks =log10_minor_break()) + theme_bw() + xlim(0,500) +
  geom_line(aes(x=x+100, y=y, colour="red")) +
  geom_line(aes(x=x+200, y=y, colour="blue")) +
  coord_cartesian(xlim=c(0,450))

print(p)
pdf("Memcurve.pdf",width=7.2,height=4)
print(p)
dev.off()
