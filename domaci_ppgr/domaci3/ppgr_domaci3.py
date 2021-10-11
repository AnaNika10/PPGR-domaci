import numpy as np
import math

def convert_to_list(matrica):
  matrica_kao_lista=[]
  for i in range(3):
    red=[]
    red.append(matrica[i,0])
    red.append(matrica[i,1])
    red.append(matrica[i,2])
    matrica_kao_lista.append(red)
  return matrica_kao_lista
def check_matrix(A):
  A=np.matrix(A)
  should_be_identity = np.allclose(A.dot(A.T), np.identity(3, np.float))
  should_be_one = np.allclose(np.linalg.det(A), 1)
  return should_be_identity and should_be_one
def norm(product):
   return np.sqrt(np.sum(product**2))

def normalize_vector(product):
  return product / norm(product)
# trazi sopstveni vektor matrice za lambda=1
def SopsVec(A):
    R=A-np.identity(3)
    product=np.cross(R[0],R[1])
    if(any(product)):
      product=normalize_vector(product)  
    return product
# fja trazi vektor koji je normalan na prosledjeni vektor
def perpendicular_vector(v):
    if v[1] == 0 and v[2] == 0:
        if v[0] == 0:
            raise ValueError('zero vector')
        else:
            return np.cross(v, [0, 1, 0])
    return np.cross(v, [1, 0, 0])

def Euler2A(fi,teta,psi):
  oko_x=np.matrix([[1,0,0],[0,math.cos(fi),-math.sin(fi)],[0,math.sin(fi),math.cos(fi)]])
  oko_y=np.matrix([[math.cos(teta),0,math.sin(teta)],[0,1,0],[-math.sin(teta),0,math.cos(teta)]])
  oko_z=np.matrix([[math.cos(psi),-math.sin(psi),0],[math.sin(psi),math.cos(psi),0],[0,0,1]])
  #A = np.dot(oko_z,np.dot(oko_y,oko_x))
  A= oko_z @ (oko_y @ oko_x)
  return A
print("Unesite Ojelerove uglove u stepenima")
fi=int(input())
teta=int(input())
psi=int(input())
fi=math.radians(fi)
teta=math.radians(teta)
psi=math.radians(psi)
print("************************Euler2A************************")
print("Matrica A")
matrica=Euler2A(fi,teta,psi)
#matrica=Euler2A(-math.atan(1/4),-math.asin(8/9),math.atan(4))
#matrica=Euler2A(math.pi/4,0,math.pi/4)
print(matrica)

def A2AxisAngle(P):
  p=SopsVec(P)
  u=perpendicular_vector(p) #normalan na p
  uprim=P @ u

  fi=math.acos(   np.sum(np.dot(u,uprim))/(norm(u)*norm(uprim))   )
  mesoviti_proizvod= np.sum(u @ (np.cross(uprim,p)))
  if(mesoviti_proizvod<0):
    p=-p
  return p,fi
#C=[[1/9,-8/9,-4/9],[4/9,4/9,-7/9],[8/9,-1/9,4/9]]
print("************************A2AxisAngle************************")
matrica_kao_lista=convert_to_list(matrica)
if check_matrix(matrica_kao_lista):
  vektor_p,ugao_fi=A2AxisAngle(matrica_kao_lista)
  
  print("Vektor p i ugao rotacije fi (u radijanima i stepenima):")
  print(vektor_p)
  print(ugao_fi,int(np.round(math.degrees(ugao_fi))))
else:
   print("Niste uneli matricu koja predstavlja kretanje")

def Rodrigez(p,fi):
  p=normalize_vector(np.array(p))
  px=np.matrix([[0,-p[2],p[1]],[p[2],0,-p[0]],[-p[1],p[0],0]])
  p=p.reshape(3,1)
  pt=p.reshape(1,3)
  ppt=np.dot(p,np.transpose(p))
  R=ppt+math.cos(fi)*(np.identity(3)-ppt)+math.sin(fi)*px
  return R
print("************************Rodrigezova f-la:************************")
print("Matrica A:")
matrica_rodrigez=Rodrigez(vektor_p,ugao_fi)
np.set_printoptions(suppress=True)
print(np.matrix(matrica_rodrigez))

def A2Euler(A):
  if A[2][0]<1 :
    if(A[2][0]>-1):
      fi=math.atan2(A[1][0],A[0][0])
      teta=math.asin(-A[2][0])
      psi=math.atan2(A[2][1],A[2][2])
    else:
      fi=math.atan2(-A[0][1],A[1][1])
      teta=math.pi/2
      psi=0
  else:
      fi=math.atan2(-A[0][1],A[1][1])
      teta=-math.pi/2
      psi=0
  return psi,teta,fi
P=[[1/9,-8/9,-4/9],[4/9,4/9,-7/9],[8/9,-1/9,4/9]]
print("************************A2Euler************************")
matrica_rodrigez_kao_lista=convert_to_list(matrica_rodrigez)
if check_matrix(matrica_rodrigez_kao_lista):
  print("Uglovi psi, teta i fi")
  psi,teta,fi=A2Euler(matrica_rodrigez_kao_lista)
  print(psi,teta,fi)
  print(int(math.degrees(psi)),int(math.degrees(teta)),int(math.degrees(fi)))
else:
   print("Niste uneli matricu koja predstavlja kretanje")

def AngleAxis2Q(p,fi):
  w=math.cos(fi/2)
  p=normalize_vector(np.array(p))
  [x,y,z]=math.sin(fi/2)*p
  q=[x,y,z,w]
  return q
print("************************AngleAxis2Q************************")
print("Kvaternion q:")
kvanterion_x,kvanterion_y,kvanterion_z,kvanterion_w=AngleAxis2Q(vektor_p,ugao_fi)
print(kvanterion_x,"i + ",kvanterion_y,"j + ",kvanterion_z,"k + ",kvanterion_w)

def Q2AxisAngle(q):
  q=normalize_vector(np.array(q))
  if(q[3]<0):
    q=-q 
  fi=2*math.acos(q[3])
  if (q[3]==1 or q[3]==-1):
    p=[1,0,0]
  else:
    p=normalize_vector(np.array([q[0],q[1],q[2]]))
  return p,fi
kvanterion_q=[kvanterion_x,kvanterion_y,kvanterion_z,kvanterion_w]
p,fi=Q2AxisAngle(kvanterion_q)
print("************************Q2AxisAngle************************")
print("Vektor p i ugao fi (u radijanama i stepenima)")
print(p)
print(fi,int(np.round(math.degrees(ugao_fi))))

