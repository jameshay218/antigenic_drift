N <- 10
rate1 <- 6
base <- ggplot(data.frame(x=c(0,N)),aes(x)) + stat_function(fun=dexp,geom="line",colour="red",args=list(rate=rate1)) +
  xlab("Size of mutation") +
  ylab("Probability (given a mutation did occur)") +
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme(
    text=element_text(colour="gray20",size=14),
    plot.title=element_text(size=28),
    legend.text=element_text(size=20,colour="gray20"),
    legend.title=element_text(size=20,colour="gray20"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line=element_line(colour="gray20"),
    axis.line.x = element_line(colour = "gray20"),
    axis.line.y=element_line(colour="gray20"),
    axis.text.x=element_text(colour="gray20"),
    panel.background=element_blank(),
    axis.text.y=element_text(colour="gray20"))
base
print(base)