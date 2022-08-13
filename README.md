<h1>Image Feature Detection</h1>

The purpose of this repository is to practically understand various algorithms and methodologies related to identifying low-level features in an image, such as corners and edges, and use that information to classify higher order image features such as lines. The techniques documented are:

<ul>
  <li> Corner detection using Harris-Stephens corner detection algorithm
  <li> Edge detection using Canny edge detection algorithm
  <li> Hough transform
  <li> Line detection and fitting
</ul>

<h2>Corner detection using Harris-Stephens corner detection algorithm</h2>

A corner, in an image, is defined as a region around which there are drastic changes in a characteristic such as intensity, color, or texture, in all directions. A popular algorithm used to detect corners in an image is the Harris-Stephens corner detection algorithm. The brilliance of this algorithm is that it takes the differential of the corner with respect to direction directly. This makes this algorithm more accurate and faster than its predecessors. The logic behind this algorithm is to consider a window around the pixel under consideration, shift this window by a small amount in all directions, and quantify the amount of change in pixel value, in each direction. The measure of change used is the Sum of Squared Difference (SSD) between the pixel values before and after the shift. For corners, those pixels are identified where the SSD is large in all 8 directions. 

![image](https://user-images.githubusercontent.com/67223688/184458903-f1c7d6e1-dd31-4d59-8b04-1c62567cb262.png)

<h2>Edge detection using Canny edge detection algorithm</h2>

While edges may be detected using the Harris-Stephen algorithm described above, a more specialized tool at disposal is the Canny edge detection algorithm. This algorithm defines an object function that is to be optimized. The algorithm uses a filter based on the first derivative of a Gaussian in order to compute the intensity of the gradients. The Gaussian reduces the effect of noise present in the image. Then, potential edges are thinned by removing non-maximum pixels of the gradient magnitude. Finally, edge pixels are kept or removed using hysteresis thresholding on the gradient magnitude. These steps can be enumerated as:

<ol>
  <li> Apply Gaussian filter to smooth the image in order to remove the noise
  <li> Find the intensity gradients of the image
  <li> Apply non-maximum suppression to get rid of spurious edges
  <li> Apply hysteresis thresholding
</ol>

![image](https://user-images.githubusercontent.com/67223688/184459038-37294376-1b84-40ca-9810-65ac7f2cb056.png)

<h2>Hough transform</h2>

The Hough transform is an algorithm developed to identify complex lines in images. While this method can be used to identify any shape that can be expressed as a model. The main intuition behind the algorithm is to parameterize the model of the shape to be identified. Doing so coverts the points of interest from an image space into a parameter space. The Hough transform takes a binary image as input. Thus, the input image has to be either the edge image attained using the Canny edge detection process described above or must be the binarized version of the original image.

![image](https://user-images.githubusercontent.com/67223688/184459098-00ad3fae-a70d-46e1-ac38-33ab41d4bc91.png)

<h2>Line detection and fitting</h2>

Once the Hough transform of an image is attained, it can be used to determine the lines present in the image.

![image](https://user-images.githubusercontent.com/67223688/184459157-f68d1f92-e7bf-4f7c-b5ac-7befcf0a29c9.png)
