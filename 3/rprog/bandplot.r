library(lattice)

bandplot <- function(formula, data, band.color = NULL, edge.color = 'black', show.whole = FALSE, ...) {
  my.panel.bands <- function(x, y, upper = data$high, lower = data$low, fill, col,
                             subscripts, ..., font, fontface)
  { 
    upper <- upper[subscripts[order(x)]]
    lower <- lower[subscripts[order(x)]]
    
    #Почему точки не отсортировали в примере?
    x <- sort(x)
    
    panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                  col = if(is.null(band.color)) fill else band.color, border = FALSE,
                  ...)
     
    panel.lines(x, upper, col = edge.color, lty = 'dashed' );
    panel.lines(x, lower, col = edge.color, lty = 'dashed');
  }
  
  xyplot(formula, data = data,
         panel = function(x, y, ...){
           panel.superpose(x, y, panel.groups = my.panel.bands, type='l', col=if(is.null(band.color)) 'gray' else band.color,...);
           # Поменял type на l с b, потому что некрасиво.
           panel.xyplot(x, y, type='l', cex=0.6, lty=1,...) 
         }
         , prepanel = if (show.whole) function(x, y, upper, lower) {
           p<-prepanel.default.xyplot(x,y,...);
           p$ylim <- c(min(lower), max(upper));
           p
         } else prepanel.default.xyplot(x,y,...)
         ,...
  ) 
}   

x = 1:10
y = x
df <- data.frame(x = x, y = y)
 
bandplot(y ~ x, df, band.color = 'green', edge.color = 'blue', show.whole = T, upper=y+2, lower=y-1, groups = (x > 0))
bandplot(y ~ x, df, edge.color = 'green', upper=y+1, lower=y-1, groups = (x < 0))
 
bandplot(y ~ x, df, upper = y^2, lower = y-1,show.whole = T, groups=(x>0))
bandplot(y ~ x, df, upper = y^2, lower = y-1, groups=(x>0))


