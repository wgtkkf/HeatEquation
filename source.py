"""
Coded by: Takuro Tokunaga
Separation of variables
Last modified: April, 19, 2017
Updated: September 7, 2026
"""
import math

class Comments:
    def __init__(self, comment1, comment2):
        self.com1 = comment1
        self.com2 = comment2

    def begin(self):
        print ("### " + str(self.com1)  + "     ###")

    def end(self):
        print ("### " + str(self.com2)  + "   ###")

class WriteToText:
    def __init__(self, filepath: str, header: str):
        self.filepath = filepath
        self.header = header

    def output(self, data_array, y_start, dy):
        with open(self.filepath, "w") as f:
            f.write(f"{self.header}\n")
            
            y = y_start
            for value in data_array:
                f.write(f"{y} {value}\n")
                y += dy

class Analytical:
    # constant parametres    
    N = 1
    NMAX = 100
    IMAX = 100
        
    # fixed temperature
    T1 = 0
    T2 = 50
    T3 = 100

    def __init__(self, w: float, l: float, x: float):
        #
        self.tsteady = [0 for i in range(Analytical.IMAX)]
        self.t1steady = [0 for i in range(Analytical.IMAX)]
        self.t2steady = [0 for i in range(Analytical.IMAX)]
        self.t3steady = [0 for i in range(Analytical.IMAX)]        

        # dimension        
        self.w = w
        self.l = l
        self.x = x
      
    def solution(self):
        # eigenvalues
        ln1 = 0
        ln2 = 0                
            
        # coordination        
        y = 0
        ymax = self.w
        dy = self.w/Analytical.IMAX
            
        # coefficient        
        c2 = (2*Analytical.T2)/math.pi                                    
        c3 = (2*Analytical.T3)/math.pi
            
        for i in range(0, Analytical.IMAX):
            for n in range(1, Analytical.NMAX):
                ln1 = (n*math.pi/self.l) # eigenvalue radian
                ln2 = (n*math.pi/self.w) # eigenvalue radian
                self.t1steady[i] = self.t1steady[i] + c3*((pow((-1),n+1)+1)/n)*math.sin(ln1*self.x)*math.sinh(ln1*y)/math.sinh(ln1*self.w)
                self.t2steady[i] = self.t2steady[i] + c3*((pow((-1),n+1)+1)/n)*math.sin(ln2*y)*math.sinh(ln2*self.x)/math.sinh(ln2*self.l)
                self.t3steady[i] = self.t3steady[i] + c2*(-math.cosh(ln2*self.l)/math.sinh(ln2*self.l))*((pow((-1),n+1)+1)/n)*math.sin(ln2*y)*\
                                        ((-math.sinh(ln2*self.l)/math.cosh(ln2*self.l))*math.cosh(ln2*self.x)+math.sinh(ln2*self.x))
                                      
            self.tsteady[i] = self.t1steady[i] + self.t2steady[i] + self.t3steady[i]            
            y = y + dy
                  
        # output for folder
        export = WriteToText("./excel/analytical.txt", "y temperature")
        export.output(self.tsteady, y_start=0, dy=dy)        

def main():
    # instance
    msg = Comments('Calculation started.', 'Calculation completed.')            
    solve = Analytical(50*pow(10,-3), 100*pow(10,-3), 0.092) # dimension of geometry, [m] and x-coordinate to be extracted [m]    

	# call methods
    msg.begin()
    solve.solution()        
    msg.end()

# --- main routine ---
if __name__ == '__main__':    
    main()