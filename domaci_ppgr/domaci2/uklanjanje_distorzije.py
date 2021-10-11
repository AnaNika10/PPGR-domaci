import cv2
import numpy as np
import math
import matplotlib.pyplot as plt
img = cv2.imread("zgrada1.jpg")
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
    return P

originals=[]
print(" Unesite koordinate na originalnoj slici koje zelite da ispravite ")
for i in range(4):
    x=float(input("Unesite x koordinatu:  "))
    y=float(input("Unesite y koordinatu:  "))
    originals.append([x,y])
originals=np.float32(originals)
images=[]
print(" Unesite koordinate provougaonika u koje preslikavate ")
for i in range(4):
    x=float(input("Unesite x koordinatu:  "))
    y=float(input("Unesite y koordinatu:  "))
    images.append([x,y])
images=np.float32(images)

#originals = np.float32([[263,949],[263,56],[663,232],[663,991]])
#images = np.float32([[263,949], [263,56], [663,56], [663,949]])

pts1 = np.array([[x, y, 1] for [x, y] in originals])
pts2 = np.array([[x, y, 1] for [x, y] in images])
M=DLT(pts1,pts2)

dst = cv2.warpPerspective(img, M, (img.shape[1],img.shape[0]))

plt.subplot(121),plt.imshow(img),plt.title('Input')
plt.subplot(122),plt.imshow(dst),plt.title('Output')
cv2.imwrite('output.jpg',dst)
plt.show()



