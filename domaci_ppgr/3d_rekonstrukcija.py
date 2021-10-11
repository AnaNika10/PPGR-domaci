import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#Fale: 8, 16, 24
x1 = [815, 110, 1]
x2 = [952, 161, 1]
x3 = [956, 125, 1]
#x3 = [989, 125, 1]
x4 = [855, 78, 1]
x5 = [790, 306, 1]
x6 = [913, 357, 1]
x7 = [950, 319, 1]
x9 = [322, 346, 1]
x10 = [453, 368, 1]
x11 = [508, 272, 1]
x12 = [385, 249, 1]
x13 = [363, 560, 1]
x14 = [477, 583, 1]
x15 = [524, 485, 1]
x17 = [136, 552, 1]
x18 = [432, 762, 1]
x19 = [818, 381, 1]
x20 = [546, 250, 1]
x21 = [175, 655, 1]
x22 = [448, 864, 1]
x23 = [804, 489, 1]
#Fale: 5, 13, 17, 21
y1 = [907, 450, 1]
y2 = [810, 560, 1]
y3 = [917, 612, 1]
y4 = [1011, 491, 1]
y6 = [774, 772, 1]
y7 = [863, 823, 1]
y8 = [957, 703, 1]
y9 = [295, 73, 1]
y10 = [251, 121, 1]
y11 = [372, 137, 1]
y12 = [413, 87, 1]
y14 = [288, 326, 1]
y15 = [396, 342, 1]
y16 = [436, 285, 1]
y18 = [143, 318, 1]
y19 = [527, 531, 1]
y20 = [741, 344, 1]
y22 = [161, 426, 1]
y23 = [532, 645, 1]
y24 = [735, 457, 1]

def getEquation(a,b):
  return np.array([a[0]*b[0],a[1]*b[0],a[2]*b[0],a[0]*b[1],a[1]*b[1],a[2]*b[1],a[0]*b[2],a[1]*b[2],a[2]*b[2] ])
jed1 = getEquation(x1,y1)
jed2 = getEquation(x2,y2)
jed3 = getEquation(x3,y3)
jed4 = getEquation(x4,y4)
jed6 = getEquation(x6,y6)
jed7 = getEquation(x7,y7)
jed9 = getEquation(x9,y9)
jed10 = getEquation(x10,y10)
jed11 = getEquation(x11,y11)
jed12 = getEquation(x12,y12)
jed14 = getEquation(x14,y14)
jed15 = getEquation(x15,y15)
jed18 = getEquation(x18,y18)
jed19 = getEquation(x19,y19)
jed20 = getEquation(x20,y20)
jed22 = getEquation(x22,y22)
jed24 = getEquation(x23,y23)
Matrica = np.array([jed1,jed2,jed3,jed4,jed6,jed7,jed9,jed10,jed11,jed12,jed14,jed15,jed18,jed19,jed20,jed22,jed24])

_, _, V = linalg.svd(Matrica)
V= V[-1:]
FF=V.reshape(3,3)
print("****F******")
print(FF)
#print("****Provera da li je det jednaka 0 (priblizno)******")
#print(np.linalg.det(FF)) 
U, D, V = linalg.svd(FF)
e1= V[-1:]
e2 = [U[0,2],U[1,2],U[2,2]]
def singularity_constraint(U, D, V):
    D1 = np.diag(np.diag([1, 1, 0]) @ D)
    F1 = U @ D1 @ V
    return F1
FF1 = singularity_constraint(U, D, V)
print("*****F1******")
print(FF1)
print("*****e1 i e2 afine koordinate*****")
e1=(e1/e1[0][2])[0]
e2=(e2/e2[2])
with np.printoptions(precision=10, suppress=True):
  print(e1)
with np.printoptions(precision=10, suppress=True):
  print(e2)
def MissingPoints(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10):
  x=np.cross(np.cross(np.cross(np.cross(a1,a2), np.cross(a3,a4)),a5), np.cross(np.cross(np.cross(a6,a7), np.cross(a8,a9)),a10))
  x = abs(x/x[2])
  x=np.round(x)
  return x
x8=MissingPoints(x5,x6,x4,x3,x7,x1,x5,x3,x7,x4)
x16=MissingPoints(x9,x10,x13,x14,x15,x13,x9,x11,x15,x12)
x24=MissingPoints(x20,x19,x21,x22,x23,x17,x21,x19,x23,x20)
y5=MissingPoints(y4,y8,y3,y7,y1,y1,y2,y8,y7,y6)  
y13=MissingPoints(y10,y14,y11,y15,y9,y10,y9,y11,y12,y14)
y17=MissingPoints(y18,y19,y22,y23,y20,y19,y20,y23,y24,y18)
y21=MissingPoints(y23,y22,y20,y17,y24,y18,y22,y20,y24,y17)
T1 = [ [1,0,0,0], [0,1,0,0], [0,0,1,0] ]
T1=np.array(T1)
print("*****T1*****")
print(T1)
def vec(p1,p2,p3):
  matr = np.array([[0,-p3,p2], [p3,0,-p1], [-p2,p1,0]])
  return matr
E2 = vec(e2[0],e2[1],e2[2])
#print("E2")
#with np.printoptions(precision=10, suppress=True):
  #print(E2)
T2=np.array(E2@FF1)
e2_novo=np.array(e2.reshape(3,1))
print("*****T2*****")
T2=np.hstack((T2,e2_novo))
with np.printoptions(precision=10, suppress=True):
  print(T2)
T2=np.array(T2)
def FourEq(x, y):
    eq = []
    eq.append(x[1]*T1[2] - x[2]*T1[1])
    eq.append(-x[0]*T1[2] + x[2]*T1[0])
    eq.append(y[1]*T2[2] - y[2]*T2[1])
    eq.append(-y[0]*T2[2] + y[2]*T2[0])
    return eq
def Affinize(x):
    x = x/x[-1]
    x = x[:-1]
    return x
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.invert_zaxis()

def Reconstruct(x, y):
    rekonstruisane = []
    for (xi, yi) in zip(x, y):
        jne = FourEq(xi, yi)
        _, _, V = np.linalg.svd(jne)
        V_nova = V[-1]
        V_nova[-1] = V_nova[-1] * 400
        rekonstruisane.append(Affinize(V_nova))
    return rekonstruisane
left = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24]
right = [y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20, y21, y22, y23, y24]
rek = Reconstruct(left, right)
def rekonstrukcija3D(rek):
    cookies =[rek[0], rek[1],rek[1],rek[2],rek[2],rek[3],rek[3],rek[0],
            rek[4],rek[5],rek[5],rek[6], rek[6],rek[7], rek[7],rek[4],
            rek[0],rek[4] ,rek[5],rek[1],rek[2],rek[6],rek[7],rek[3]]

    tea =[ rek[8], rek[9],rek[9],rek[10],rek[10],rek[11],rek[11],rek[8],
            rek[12],rek[13],rek[13],rek[14], rek[14],rek[15], rek[15],rek[12],
            rek[8],rek[12],rek[13],rek[9],rek[10],rek[14],rek[15],rek[11]
            ]
           
    router =[ rek[16], rek[17],rek[17],rek[18],rek[18],rek[19],rek[19],rek[16],
              rek[20],rek[21],rek[21],rek[22], rek[22],rek[23], rek[23],rek[20],
              rek[16],rek[20],rek[21],rek[17],rek[18],rek[22],rek[23],rek[19] ]       
    ax.view_init(20,5)
    for (k, boja) in zip([cookies, tea, router], ['red', 'green', 'blue']):
        x=[]
        y=[]
        z=[]
        for e in k:
          x.append(e[0])
        for e in k:
          y.append(e[1])
        for e in k:
          z.append(e[2])
        ax.plot3D(x, y, z, boja)
    plt.show()
rekonstrukcija3D(rek)