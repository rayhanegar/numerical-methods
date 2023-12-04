er = lambda cc, cp: abs((cc-cp)/cc) * 100 

# Bisection method
def bisection(a, b, threshold, fx):

    # calculate initial c0
    c0 = (a+b) / 2

    print("Init: ", a, b, c0)
    print("a\tb\tfx(a)\tfx(b)\tc0\tfx(c0)\tc1\tfx(c1)\ter")

    condition0 = fx(a) * fx(b)
    
    # first iteration outside loop to find c1
    if (condition0 < 0):
        if ((fx(a) * fx(c0)) > 0):
            a = c0
        elif ((fx(a) * fx(c0)) < 0):
            b = c0
        else:
            return b
    
    # calculate c1
    c1 = (a+b) / 2
    
    print(f'{a:.3f}\t{b:.3f}\t{fx(a):.3f}\t{fx(b):.3f}\t{c0:.3f}\t{fx(c0):.3f}\t{c1:.3f}\t{fx(c1):.3f}\t{er(c0, c1):.3f}')

    while(er(c0, c1) > threshold):
        
        # calculate current condition
        condition_current = fx(a) * fx(c1)

        # update either a or b with c1
        if(condition_current > 0):
            a = c1
        elif(condition_current < 0):
            b = c1

        # update c0 with c1
        c0 = c1

        # calculate new c1 using recently updated a and b
        c1 = (a+b) / 2

        print(f'{a:.3f}\t{b:.3f}\t{fx(a):.3f}\t{fx(b):.3f}\t{c0:.3f}\t{fx(c0):.3f}\t{c1:.3f}\t{fx(c1):.3f}\t{er(c0, c1):.3f}')
    
    return c1

#Regula Falsi method
def regula_falsi(a, b, threshold, fx):

    print("Init: ", a, b)
    print("a\tb\tfx(a)\tfx(b)\tc0\tfx(c0)\tc1\tfx(c1)\ter abs(fx(c0))")

    if fx(a) * fx(b) < 0:

        # calculate c0
        c0 = b - (fx(b) * (b - a) / (fx(b) - fx(a)))

        while abs(fx(c0)) > threshold:
            
            # calculate c1, first iteration result is identical with first c0
            c1 = b - (fx(b) * (b - a)/(fx(b) - fx(a)))

            # update a or b with c1
            if fx(a) * fx(c1) > 0:
                a = c1
            elif fx(a) * fx(c1) < 0:
                b = c1
            
            print(f'{a:.3f}\t{b:.3f}\t{fx(a):.3f}\t{fx(b):.3f}\t{c0:.3f}\t{fx(c0):.3f}\t{c1:.3f}\t{fx(c1):.3f}\t{abs(fx(c0)):.3f}')

            # update c0 with c1
            c0 = c1
            
        return c1

#Newton-Rhapson method
def newton_raphson(init_guess, fx, fprimex, max_iteration, threshold):
    
    next = lambda x0, fx, fprimex: x0 - (fx(x0) / fprimex(x0))
    x0 = init_guess

    print("Init: ", x0, "Max-Iter: ", max_iteration, "Threshold: ", threshold)
    print("i\tx0\tfx(0)\tf'x(x0)\tx1\ter")

    for i in range(max_iteration):

        # calculate value for initial guess
        fx0 = fx(x0)

        # calculate value for f'(x0)
        fprimex0 = fprimex(x0)

        # calculate next point, x1
        x1 = next(x0, fx, fprimex)

        # return x1 if less than threshold
        if abs(x0 - x1) < threshold:
            return x1
        
        print(f'{i}\t{x0:.3f}\t{fx0:.3f}\t{fprimex0:.3f}\t{x1:.3f}\t{abs(x0-x1):.3f}')
        
        # update x0 to x1
        x0 = x1

#Secant method
def secant(a, b, threshold, fx, max_iter=20):
    
    diff = lambda a, b, fx: fx(b) * (b-a) / (fx(b) - fx(a))
    next = lambda a, b, fx: b - (fx(b) * (b-a) / (fx(b) - fx(a))) 

    i = 1
    er = abs(diff(a, b, fx))

    print('i\ta\tb\tfx(a)\tfx(b)\ter')

    while (er >= threshold) & (i <= max_iter):

        print(f'{i:d}\t{a:.3f}\t{b:.3f}\t{fx(a):.3f}\t{fx(b):.3f}\t{er:.3f}')
        
        # calculate next b value
        next_b = next(a, b, fx)

        # update a and b values
        a = b
        b = next_b

        # calculate error using recently updated a and b
        er = abs(diff(a, b, fx))

        i+=1
    
    return b