'''
P11=[595,301,1]
P21=[292,517,1]
P31=[157,379,1]
P51=[665,116,1]
P61=[304,295,1]
P71=[135,163,1]
P81=[509,43,1]
'''
import numpy as np

def Afinize(P):
    P_afino=[]
    P_afino.append(np.round(P[0]/P[2]))
    P_afino.append(np.round(P[1]/P[2]))
    P_afino.append(1)
    return P_afino

def Nevidljivo(P1,P2,P3,P5,P6,P7,P8):

    Xb1=np.round(Afinize(np.cross(np.cross(P2,P6),np.cross(P1,P5))))
    Xb2=np.round(Afinize(np.cross(np.cross(P2,P6),np.cross(P3,P7))))
    Xb3=np.round(Afinize(np.cross(np.cross(P1,P5),np.cross(P3,P7))))

    XbZbir=Xb1+Xb2+Xb3
    XbTeziste=[x/3 for x in XbZbir]
    Xb=Afinize(np.round(XbTeziste))
    print(Xb) #[252.0, 2247.0, 1]
    Yb1=np.round(Afinize(np.cross(np.cross(P5,P6),np.cross(P7,P8))))
    Yb2=np.round(Afinize(np.cross(np.cross(P2,P1),np.cross(P7,P8))))
    Yb3=np.round(Afinize(np.cross(np.cross(P2,P1),np.cross(P5,P6))))
    YbZbir=Yb1+Yb2+Yb3
    YbTeziste=[y/3 for y in YbZbir]
    Yb=Afinize(np.round(YbTeziste))
    print(Yb)#[2221.0, -106.0, 1]
    P4=np.cross(np.cross(P8,Xb),np.cross(P3,Yb))
    print('Homogene')
    print(P4) #[-1.95524056e+09 -2.34987757e+09 -4.31770800e+06]
    P4=Afinize(P4)
    P4_afino=Afinize(P4)
    return P4_afino

P1=[646,657,1]
P2=[295,840,1]
P3=[100,674,1]
P5=[763,213,1]
P6=[308,309,1]
P7=[56,181,1]
P8=[503,119,1]
P4_afine=Nevidljivo(P1,P2,P3,P5,P6,P7,P8)
print('Afinne')
print(P4_afine)#[453.0, 544.0, 1]  Trazena tacka je (453,544)

