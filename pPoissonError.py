from pPoisson import *
wStep = 0.01
for iNum in range(1, 20):
    print "-----------------------------"
    print "Poissonian error of", iNum, "event(s):"
    endLowErr = 0
    step = iNum-wStep/2.
    cumLow = 0
    while endLowErr==0 and step>0:
        fc = fcPoisson(step)
        cumLow = cumLow + fc.getTerm(iNum)*wStep
        if cumLow >= 0.34:
            endLowErr = step
        else:
            step = step-wStep
    endUpErr = -1
    step = iNum+wStep/2.
    cumUp = 0
    while endUpErr==-1:
        fc = fcPoisson(step)
        cumUp = cumUp+ fc.getTerm(iNum)*wStep
        if cumUp >= 0.34:
            endUpErr = step
        else:
            step = step+wStep
    print iNum, "+", endUpErr-iNum, "-", iNum-endLowErr
    print "cf.", iNum, "+", 0.5 + math.sqrt(iNum+0.25), "-", -0.5 + math.sqrt(iNum+0.25)
