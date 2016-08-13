import math

class fcPoisson:
    def __init__(self, lam):
        self.lam = lam
    def getParameter(self):
        return self.lam
    def getTerm(self, iTerm):
        return math.pow(self.lam, iTerm)*math.exp(-self.lam)/float(math.factorial(iTerm))
    def getCum(self, iStart, iEnd):
        sumTerm = 0
        for jTerm in range(iStart, iEnd):
            sumTerm = sumTerm + self.getTerm(jTerm)
        return sumTerm
    def getCumZero(self, iEnd):
        return self.getCum(0, iEnd)
