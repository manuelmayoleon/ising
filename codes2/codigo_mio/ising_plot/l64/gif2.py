import imageio
import os
from pathlib import Path
# You must create a folder named "source_images" and put the pictures inside 
image_path = Path('source_images')
images = list(image_path.glob('*.png'))
#sorting files to put in correct order to convert to .gif 
images = sorted(images)

image_list = []

# print(images)
for file_name in images:
    image_list.append(imageio.imread(file_name))

imageio.mimwrite('animated_from_images.gif', image_list,fps=10)
imageio.mimread('animated_from_images.gif')
