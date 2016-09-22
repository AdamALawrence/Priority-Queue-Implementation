DISTANCE FIELDS & FAST MARCHING METHOD


NAME:  Adam Lawrence


COLLABORATORS AND OTHER RESOURCES:
Salles
Toshi
Aritra C.


ESTIMATE OF # OF HOURS SPENT ON THIS: 9



NAIVE ALGORITHM:

Order Notation: O( i*j*i2*j2 )

Timing Experiment Data:
squiggle_30x30.ppm out.ppm naive_method greyscale
real	0m0.013s
user	0m0.010s
sys	0m0.002s

data_structures_300x300.ppm out.ppm naive_method rainbow
real	0m11.902s
user	0m11.878s
sys	0m0.018s

Discussion:
This algorithm, although simple, has a high order notation because it compares every pixel to every other pixel. In a 30x30 image thats 1800 checks. The larger the image, the run time exponentially increases



IMPROVED ALGORITHM:

k = # of black pixels
Order Notation: O( (i*j) + (i*j*k) ) => O( i*j*k )

Timing Experiment Data:
squiggle_30x30.ppm out.ppm improved_method rainbow
real	0m0.005s
user	0m0.002s
sys	0m0.001s

data_structures_300x300.ppm out.ppm improved_method rainbow
real	0m8.317s
user	0m8.304s
sys	0m0.009s

Discussion:
The improved method first scans for all the black pixels in the image
This eliminates the need to compare all the pixels to each other.
After storing all the black pixels, I then compare each pixel to each black pixel and find the max distance that way. This results in a triple loop rather than a quad loop, improving the run time.



FAST MARCHING METHOD (with a map):

k = # of black pixels
v = size of values
h = average size of heap
Order Notation: O( i*j[find black pxls] + klogv[add black pixels to queue] + hlogh[push and update functions] + logh[pop from queue] )

Timing Experiment Data:
squiggle_30x30.ppm out.ppm pq_with_map rainbow
real	0m0.004s
user	0m0.002s
sys	0m0.001s

data_structures_300x300.ppm out.ppm pq_with_map rainbow
real	0m0.172s
user	0m0.166s
sys	0m0.003s

Discussion:
Even though the Fast Marching Method contains more code, it is exponentially faster than the other two methods. The functions that use percolate up or down have a log run time.
This results in chaining of logs, each with its own purpose, whether it be updating, popping, or finding pixels within the queue/heap.


DISTANCE FIELD VISUALIZATIONS FOR EXTRA CREDIT:
I changed the way colors were represented.
I changed colors that would display a pleasant holiday feeling from a foreign world - a joyous extraterrestrial christmasy feeling.
That is all.

(Image is in zip)



