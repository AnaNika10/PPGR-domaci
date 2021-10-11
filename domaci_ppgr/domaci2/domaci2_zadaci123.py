import numpy as np
def SolveSystem(a,b,c,d):
    A = np.array([    [a[0],b[0],c[0]]   , [a[1],b[1],c[1]] , [a[2],b[2],c[2]]   ])
    B = np.array([  d[0],d[1],d[2]  ])
    alfa,beta,gama = np.linalg.solve(A,B)
    return alfa,beta,gama
def Naive(a,b,c,d,ap,bp,cp,dp):
    alfa,beta,gama=SolveSystem(a,b,c,d)
    
    p1_prva_kolona=np.array([alfa*x for x in a])
    p1_druga_kolona=np.array([beta*x for x in b])
    p1_treca_kolona=np.array([gama*x for x in c])

    alfaprim,betaprim,gamaprim=SolveSystem(ap,bp,cp,dp)
    
    p2_prva_kolona=np.array([alfaprim*x for x in ap])
    p2_druga_kolona=np.array([betaprim*x for x in bp])
    p2_treca_kolona=np.array([gamaprim*x for x in cp])

    p1=np.column_stack((p1_prva_kolona,p1_druga_kolona,p1_treca_kolona))
    p2=np.column_stack((p2_prva_kolona,p2_druga_kolona,p2_treca_kolona))
    P=np.matrix.round(p2.dot(np.linalg.inv(p1)),8) 
    return P

def TwoEquations(a,ap):
    v1=np.array([0,0,0,-ap[2]*a[0],-ap[2]*a[1],-ap[2]*a[2],ap[1]*a[0],ap[1]*a[1],ap[1]*a[2]])
    v2=np.array([ap[2]*a[0],ap[2]*a[1],ap[2]*a[2],0,0,0,-ap[0]*a[0],-ap[0]*a[1],-ap[0]*a[2]])
    return v1,v2
def DLT(points_original,points_images):
    n=len(points_original)
    e1,e2=TwoEquations(points_original[0],points_images[0])
    A = np.vstack((e1,e2))
    for i in range(1,n):
        e1,e2=TwoEquations(points_original[i],points_images[i])
        A = np.vstack([A, e1])
        A = np.vstack([A, e2])
    _, _, V = np.linalg.svd(A)
    P = V[-1, :]
    P=P.reshape(3,3)
    P=np.matrix.round(P,5)
    return P
def DLTModification(points_original,points_images):
    points_original_normalized,T=Normalize(points_original)
    points_images_normalized,Tp=Normalize(points_images)
    Pdlt_normalized=DLT(points_original_normalized,points_images_normalized)
    ### DLT SA NORMALIZACIJOM KAD VRATIMO U POCETNI KOORDINATNI SISTEM #####
    TpI=np.linalg.inv(Tp)
    Pdlt_normalized_withoutT=np.matrix.round(TpI.dot(Pdlt_normalized).dot(T),5) 
    return Pdlt_normalized_withoutT
def Afinize(P):
    P_afino=[]
    P_afino.append(P[0]/P[2])
    P_afino.append(P[1]/P[2])
    P_afino.append(1)
    P_afino=np.array(P_afino)
    return P_afino
def Transformation(C1,C2,r):
    G=np.identity(3)
    G[0,2]=-C1
    G[1,2]=-C2
    S=np.identity(3)
    S[0,0]=np.sqrt(2)/r
    S[1,1]=np.sqrt(2)/r
    T=np.matmul(S,G)
    return T
def Normalize(points):
    n=len(points)
    afinized_points=[Afinize(x) for x in points]
    C1=(np.sum([x[0] for x in afinized_points]))/n
    C2=(np.sum([x[1] for x in afinized_points]))/n
    first_coord=[x[0]-C1 for x in afinized_points]
    second_coord=[x[1]-C2 for x in afinized_points]
    distances_from_O=[np.sqrt(x**2+y**2) for x,y in zip(first_coord,second_coord)]
    r=np.sum(distances_from_O)/n
    T=Transformation(C1,C2,r)
    points_normalized=[np.matrix.round(np.matmul(T,x),6) for x in points]
    points_normalized=tuple(points_normalized)
    return points_normalized,T
'''
p_given=[[0,3,5],[4,0,0],[-1,-1,6]]
a = (-3, 2, 1)
b = (-2, 5, 2)
c = (1, 0, 3)
d = (-7, 3, 1)
e = (2, 1, 2)
f = (-1, 2, 1)
g = (1, 1, 1)

ap=tuple(np.dot(p_given,a))
bp=tuple(np.dot(p_given,b))
cp=tuple(np.dot(p_given,c))
dp=tuple(np.dot(p_given,d))
ep=tuple(np.dot(p_given,e))
fp=tuple(np.dot(p_given,g))
gp=(8.02,4,4)
'''
'''
a = (1, 1, 1)
b = (5, 2, 1)
c = (6, 4, 1)
d = (-1, 7, 1)
ap = (0, 0, 1)
bp = (10, 0, 1)
cp = (10, 5, 1)
dp = (0, 5, 1)
e = (3 , 1, 1)
'''

a= (-3, -1, 1)
b = (3, -1, 1)
c = (1, 1, 1)
d = (-1, 1, 1)
ap = (-2, -1, 1)
bp = (2, -1, 1)
cp = (2, 1, 1)
dp = (-2, 1, 1)
e = (1, 2, 3)
f = (-8, -2, 1)

##################### NAIVNI ##########################
P=Naive(a,b,c,d,ap,bp,cp,dp)
print("P NAIVNI ALGORITAM")
print(P)
##################### DODAJEMO TACKE ####################
points_original=[]
points_images=[]
points_original.extend([a,b,c,d])
points_images.extend([ap,bp,cp,dp])
#ep=tuple(np.matrix.round(Afinize(np.dot(P,e))))
#fp=tuple(np.matrix.round(Afinize(np.dot(P,f))))
ep=tuple(np.dot(P,e))
fp=tuple(np.dot(P,f))
##################### DLT ############################
Pdlt4=DLT(points_original,points_images)
print("P DLT ALGORITAM SA 4 TACKE")
print(Pdlt4)
points_original.append(e)
points_images.append(ep)
points_original.append(f)
points_images.append(fp)
#points_original.append(g)
#points_images.append(gp)
Pdlt=DLT(points_original,points_images)
print("P DLT ALGORITAM SA VISE TACAKA")
print(Pdlt)
##################### POREDJENJE DLT I NAIVNI ##########################
prvi_elem_dlt=Pdlt4[2,2]
prvi_elem_naive=P[2,2]
poredjenje=Pdlt4/prvi_elem_dlt*prvi_elem_naive
print("POREDJENJE DLT I NAIVNI")
print(poredjenje)
##################### MODIFIKOVANI DLT (SA NORMALIZACIJOM)  ##########################
Pdlt_normalized_withoutT=DLTModification(points_original,points_images)
print("Pdlt_normalized (without T and Tp)")  
print(Pdlt_normalized_withoutT)
#################### POREDJENJE DLT I DLT SA NORMALIZACIJOM ##########################
prvi_elem_dlt_mofication=Pdlt_normalized_withoutT[2,2]
prvi_elem_dlt=Pdlt[2,2]
poredjenje=np.matrix.round(Pdlt_normalized_withoutT/prvi_elem_dlt_mofication*prvi_elem_dlt,6)
print("POREDJENJE DLT I DLT SA NORMALIZACIJOM")
print(poredjenje)
############################### PERMUTOVANE TACKE DLT ALGORITAM ########################
points_original_permutation=[]
points_images_permutation=[]
points_original_permutation.extend([a,c,b,d,e])
points_images_permutation.extend([ap,cp,bp,dp,ep])
Pdlt_permutaion=DLT(points_original_permutation,points_images_permutation)
print("P DLT ALGORITAM SA 5 TACAKA KAD SU PERMUTOVANE 2 I 3 (ZAKLJUCAK-NIJE OSETLJIV NA OVO)")
print(Pdlt_permutaion)
############################### TRANSFORMISANE TACKE DLT ALGORITAM ########################
C1 = [[0, 1, 2], [-1, 0, 3], [0, 0, 1]]
C2 = [[1, -1, 5], [1, 1, -2], [0, 0, 1]]
points_original_transformed=[np.dot(C1,x) for x in points_original]
points_images_transformed=[np.dot(C2,x) for x in points_images]
Pdlt_transformed=DLT(points_original_transformed,points_images_transformed)
print("P DLT ALGORITAM SA 5 TACAKA KAD SU TRANSFORMISANE")
print(Pdlt_transformed)
Pdlt_old=np.matrix.round(np.linalg.inv(C2).dot(Pdlt_transformed).dot(C1),8)
print("Proveravamo da li je rezultat DLP algoritma primenjenog na nove koordinate, isti kao rezultat starog u novim koordinatama")
print(Pdlt_old)
Pdlt_old_scaled=np.matrix.round(Pdlt_old*(Pdlt[0,0]/Pdlt_old[0,0]),7)
print("Posle skaliranja (ZAKLJUCAK-OSETLJIV JE NA OVO)")
print(Pdlt_old_scaled)
############################### TRANSFORMISANE TACKE MODIFIKOVAN DLT ALGORITAM ########################
Pdlt_transformed_normalized=DLTModification(points_original_transformed,points_images_transformed)
print("P DLT MODIFIKOVAN ALGORITAM SA 5 TACAKA KAD SU TRANSFORMISANE")
print(Pdlt_transformed_normalized)
Pdlt_old_normalized=np.matrix.round(np.linalg.inv(C2).dot(Pdlt_transformed_normalized).dot(C1),8)
print("Proveravamo da li je rezultat DLP algoritma primenjenog na nove koordinate, isti kao rezultat starog u novim koordinatama")
print(Pdlt_old_normalized)
Pdlt_old_normalized_scaled=np.matrix.round(Pdlt_old_normalized*(Pdlt_normalized_withoutT[0,0]/Pdlt_old_normalized[0,0]),5)
print("Posle skaliranja (ZAKLJUCAK- NIJE OSETLJIV JE NA OVO)")
print(Pdlt_old_normalized_scaled)

