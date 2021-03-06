#pip install opencv-contrib-python==3.4.2.16
import cv2
import numpy as np
import matplotlib.pyplot as plt

img_ = cv2.imread("2.jpeg")
img1 = cv2.cvtColor(img_,cv2.COLOR_BGR2GRAY)
img = cv2.imread("1.jpeg")
img2 = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
sift = cv2.xfeatures2d.SIFT_create()
# uzimamo zajednicke tacke na obe slike
kp1, des1 = sift.detectAndCompute(img1,None)
kp2, des2 = sift.detectAndCompute(img2,None)
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1,des2, k=2)

good = []
for m in matches:
    if m[0].distance < 0.5*m[1].distance:
        good.append(m)
matches = np.array(good)
if len(matches[:,0]) >= 4:
  org_points = np.float32([ kp1[m.queryIdx].pt for m in matches[:,0] ]).reshape(-1,1,2)
  img_points = np.float32([ kp2[m.trainIdx].pt for m in matches[:,0] ]).reshape(-1,1,2)
  ProjectionMatrix, _ = cv2.findHomography(org_points, img_points, cv2.RANSAC, 5.0)
else:
  print("Slike nemaju dovoljno zajednickih tacaka")

# uklanja crni okvir oko slike ukoliko nastane
def trim(frame):
    #crop top
    if not np.sum(frame[0]):
        return trim(frame[1:])
    #crop bottom
    elif not np.sum(frame[-1]):
        return trim(frame[:-2])
    #crop left
    elif not np.sum(frame[:,0]):
        return trim(frame[:,1:]) 
    #crop right
    elif not np.sum(frame[:,-1]):
        return trim(frame[:,:-2])    
    return frame

dst = cv2.warpPerspective(img_,ProjectionMatrix,(img.shape[1] + img_.shape[1], img.shape[0]))
dst[0:img.shape[0], 0:img.shape[1]] = img
cv2.imwrite('output.jpg',trim(dst))
plt.imshow(trim(dst))
plt.show()

