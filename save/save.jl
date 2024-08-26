# ----------------------------------------------------
# file zermelo loss1

# second choise of regularization 
#function indicator(x, a, b)
#    k = 120 
#    g1 = 1 / (1 + exp(-k * (x - a)))  
#    g2 = 1 / (1 + exp(-k * (b - x)))  
#    return g1 * g2
#end
#fNC(x)  = indicator(x, 0.5, 3.5) ;
#fC(x)  = (1. - fNC(x));

# ----------------------------------------------------
# file zermelo loss2

#function indicator(x, a, b)
#    k = 17
#    g1 = 1 / (1 + exp(-k * (x - a)))  
#    g2 = 1 / (1 + exp(-k * (b - x)))  
#    return g1 * g2
#end
#f1(x)  = indicator(x,  4.8, 10.2) ;
#f2(x)  = indicator(x, 19.8, 25.2) ;
#fC1(x)  = (1. - f1(x))*(1. - f2(x))
#fNC1(x) = (f2(x) + f1(x))
#plot(fNC1, -1, 30);
