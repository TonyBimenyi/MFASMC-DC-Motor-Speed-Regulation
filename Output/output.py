import numpy as np

class Output:
    def __init__(self, rT = 1040, m=2)
        self.rT = rT
        self.m = m

    def update_y(self,k):
        if k == 1:
            return 0.51, 2.5, 3.5, 4
        else:
            y1_next = self.m/(self.rT * 2)
            y2_next = self.m/(self.rT * 3)
            y3_next = self.m/(self.rT * 4)
            y4_next = self.m/(self.rT * 5)

            return y1_next,y2_next,y3_next,y4_next