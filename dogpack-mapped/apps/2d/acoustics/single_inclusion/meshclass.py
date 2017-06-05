class cell:
   "store point values and coefficients"
   x=[]
   y=[]
   coeff=[]
   order=0
   def __init__(self,order,x,y,q):
      self.order=order
      self.x=x
      self.y=y
      self.coeff=q
     
class meshobject:
    "store cells" 
    cells={}
    def setcell(self,i,j,order,x,y,q):
       c1=cell(order,x,y,q)
       self.cells[(i,j)]=c1
