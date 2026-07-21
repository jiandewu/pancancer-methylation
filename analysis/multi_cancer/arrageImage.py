  
library(reticulate)
use_python("/anaconda3/bin/python")

# code for displaying multiple images in one figure
  
#import libraries
import cv2
from matplotlib import pyplot as plt
  
# create figure
fig = plt.figure(figsize=(10, 7))
  
# setting values to rows and column variables
rows = 2
columns = 2
  
# reading images
Image1 = cv2.imread('./multi_cancer/BRCA/figure/Tumor.vs.Ctrl_DMP_thr_p1e-15_dB0.3_heatmap.png')
Image2 = cv2.imread('./multi_cancer/COAD/figure/Tumor.vs.Ctrl_DMP_thr_p1e-15_dB0.3_heatmap.png')
Image3 = cv2.imread('./multi_cancer/LUAD/figure/Tumor.vs.Ctrl_DMP_thr_p1e-15_dB0.3_heatmap.png')
Image4 = cv2.imread('./multi_cancer/PRAD/figure/Tumor.vs.Ctrl_DMP_thr_p1e-15_dB0.3_heatmap.png')

# Adds a subplot at the 1st position
fig.add_subplot(rows, columns, 1)
  
# showing image
plt.imshow(Image1)
plt.axis('off')
plt.title("First")
  
# Adds a subplot at the 2nd position
fig.add_subplot(rows, columns, 2)
  
# showing image
plt.imshow(Image2)
plt.axis('off')
plt.title("Second")
  
# Adds a subplot at the 3rd position
fig.add_subplot(rows, columns, 3)
  
# showing image
plt.imshow(Image3)
plt.axis('off')
plt.title("Third")
  
# Adds a subplot at the 4th position
fig.add_subplot(rows, columns, 4)
  
# showing image
plt.imshow(Image4)
plt.axis('off')
plt.title("Fourth")

plt.show()
