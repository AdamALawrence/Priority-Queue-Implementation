#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>

#include "image.h"
#include "priority_queue.h"

// ===================================================================================================

// distance field method functions
double NaiveDistanceFieldMethod(Image<Color> &input, Image<DistancePixel> &distance_image);
double ImprovedDistanceFieldMethod(Image<Color> &input, Image<DistancePixel> &distance_image);
double FastMarchingMethod(Image<Color> &input, Image<DistancePixel> &distance_image);

// visualization style helper functions
Color Rainbow(double distance, double max_distance);
Color GreyBands(double distance, double max_distance, int num_bands);
Color AlienChristmas(double distance, double max_value, int num_bands);
// ===================================================================================================

int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " input.ppm output.ppm distance_field_method visualization_style" << std::endl;
    exit(1);
  }

  // open the input image
  Image<Color> input;
  if (!input.Load(argv[1])) {
    std::cerr << "ERROR: Cannot open input file: " << argv[1] << std::endl;
    exit(1);
  }

  // a place to write the distance values
  Image<DistancePixel> distance_image;
  distance_image.Allocate(input.Width(),input.Height());

  // calculate the distance field (each function returns the maximum distance value)
  double max_distance = 0;
  if (std::string(argv[3]) == std::string("naive_method")) {
    max_distance = NaiveDistanceFieldMethod(input,distance_image);
  } else if (std::string(argv[3]) == std::string("improved_method")) {
    max_distance = ImprovedDistanceFieldMethod(input,distance_image);
  } else if (std::string(argv[3]) == std::string("pq_with_map")) {
    max_distance = FastMarchingMethod(input,distance_image);
  } else if (std::string(argv[3]) == std::string("pq_with_hash_table")) {
    // EXTRA CREDIT: implement FastMarchingMethod with a hash table
  } else {
    std::cerr << "ERROR: Unknown distance field method: " << argv[3] << std::endl;
    exit(1);
  }

  // convert distance values to a visualization
  Image<Color> output;
  output.Allocate(input.Width(),input.Height());
  for (int i = 0; i < input.Width(); i++) {
    for (int j = 0; j < input.Height(); j++) {
      double v = distance_image.GetPixel(i,j).getValue();
      if (std::string(argv[4]) == std::string("greyscale")) {
	output.SetPixel(i,j,GreyBands(v,max_distance*1.01,1));
      } else if (std::string(argv[4]) == std::string("grey_bands")) {
	output.SetPixel(i,j,GreyBands(v,max_distance,4));
      } else if (std::string(argv[4]) == std::string("rainbow")) {
	output.SetPixel(i,j,Rainbow(v,max_distance));
      } else {
	// EXTRA CREDIT: create other visualizations 
          output.SetPixel(i,j,AlienChristmas(v,max_distance,4));
      }
    }
  }
  // save output
  if (!output.Save(argv[2])) {
    std::cerr << "ERROR: Cannot save to output file: " << argv[2] << std::endl;
    exit(1);
  }

  return 0;
}

// ===================================================================================================
// HELPER FUNCTIONS

double findDistance(DistancePixel* pix1, DistancePixel* pix2) {
    
    double distance;
    
    double pix1_x = pix1->getX(); double pix1_y = pix1->getY();
    double pix2_x = pix2->getX(); double pix2_y = pix2->getY();
    
    distance = sqrt((pix2_x-pix1_x)*(pix2_x-pix1_x) + (pix2_y-pix1_y)*(pix2_y-pix1_y));
    
    return distance;
}

// THIS FUNCTION CHECKS THE VALUE OF THE NEIGHBORING PIXEL SHOULD BE UPDATED
bool ShouldBeUpdated(DistancePixel* current, DistancePixel* neighbor) {
    
    double c_D = current->getValue();
    double n_D = neighbor->getValue();
    
    double dis = findDistance(current, neighbor);
    
    if ( (c_D + dis) < n_D) {
        return true;
    }
    else
        return false;
    
}

// ===================================================================================================
double NaiveDistanceFieldMethod(Image<Color> &input, Image<DistancePixel> &distance_image) {
    int w = input.Width();
    int h = input.Height();
    // return the maximum distance value
    double answer = 0;
    // loop over the pixels in the input image
    for (int i = 0; i < w; i++)  {
        for (int j = 0; j < h; j++) {
            double closest = -1;
            // loop over all other pixels in the input image
            for (int i2 = 0; i2 < w; i2++)  {
                for (int j2 = 0; j2 < h; j2++) {
                    const Color& c = input.GetPixel(i2,j2);
                    // skip all pixels that are not black
                    if (!c.isBlack()) continue;
                    // calculate the distance between the two pixels
                    double distance = sqrt((i-i2)*(i-i2) + (j-j2)*(j-j2));
                    // store the closest distance to a black pixel
                    if (closest < 0 || distance < closest) {
                        closest = distance;
                    }
                }
            }
            assert (closest >= 0);
            answer = std::max(answer,closest);
            // save the data to the distance image
            DistancePixel& p = distance_image.GetPixel(i,j);
            p.setValue(closest);
        }
    }
    return answer;
}

double ImprovedDistanceFieldMethod(Image<Color> &input, Image<DistancePixel> &distance_image) {

    // instead of comparing every pixel to every other pixel, first you find all the black pixels
    //    then you compair all other pixels to each black pixel, finding the smallest distance as you go
    double answer = 0;
    int w = input.Width();
    int h = input.Height();
    std::vector<DistancePixel> vecBlackPixels;
    
    // get the locations off all the black pixels
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            const Color& c = input.GetPixel(i, j);
            if (c.isBlack()) {
                DistancePixel p;
                p.setX(i);
                p.setY(j);
                vecBlackPixels.push_back(p);
            }
        }
    }
    // find the shortest distance between a black pixel and a non-black pixel
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
           // if (c.isBlack()) {continue;}
            //double distance = 0;
            double closest = -1;
            for (int k = 0; k < vecBlackPixels.size(); k++) {
                
                double x = vecBlackPixels[k].getX();
                double y = vecBlackPixels[k].getY();
                double distance = sqrt((x-i)*(x-i) + (y-j)*(y-j));
                
                if (closest < 0 || distance < closest) {
                    closest = distance;
                    assert(closest >= 0);
                    DistancePixel& p = distance_image.GetPixel(i, j);
                    p.setValue(closest);
                }
            }
            answer = std::max(answer, closest);
        }
    }
    
    
    return answer;
}

double FastMarchingMethod(Image<Color> &input, Image<DistancePixel> &distance_image) {
    
    int w = input.Width();
    int h = input.Height();
    
    std::vector<DistancePixel*> vecBlackPixels;
    
    for (int i = 0; i < w; i++) {                                           // get the locations off all the black pixels
        for (int j = 0; j < h; j++) {
            const Color& c = input.GetPixel(i, j);
            if (c.isBlack()) {
                DistancePixel* temp = &distance_image.GetPixel(i, j);       // set black pixels to i,j,0
                temp->setX(i);
                temp->setY(j);
                temp->setValue(0);
                vecBlackPixels.push_back(temp);
            }
            else {
                DistancePixel* temp = &distance_image.GetPixel(i, j);       // set nonblacks to i,j,reallyyybignumber
                temp->setX(i);
                temp->setY(j);
            }
        }
    }
    DistancePixel_PriorityQueue priorityQueue(vecBlackPixels);
                                                                            // intialize beginning priority queue vector with
                                                                            // only black pixels
    
    // ===================================================================================================
    // MAIN LOOP THAT GOES THROUGH QUEUE UPDATES BACKPOINTERS AND
    
    const double sq2 = sqrt(2);
    double max_distance = 0;
    while (!priorityQueue.empty()) {

        DistancePixel* currentPixel = priorityQueue.top();
        
        int c_x = currentPixel->getX();                                         // 'c_x' -> current's x coordinate
        int c_y = currentPixel->getY();                                         // 'c_y' -> current's y coordinate
        
        
        // CHECK BOUNDARIES STARTING AT LEFT NEIGHBOR GOING CLOCKWISE
        if (c_y-1 >= 0) {                                                       // check left
            
            DistancePixel* neighbor = &distance_image.GetPixel(c_x, c_y-1);     // create pixel to look at
            
            if (ShouldBeUpdated(currentPixel, neighbor)) {                      // check if pixel should be updated
                
                neighbor->setX(currentPixel->getX());
                neighbor->setY(currentPixel->getY()-1);
                neighbor->setValue(currentPixel->getValue() + 1);               // update pixel's distance value
                
                if (!priorityQueue.in_heap(neighbor)) {
                    priorityQueue.push(neighbor);                               // if pixel isn't in queue, add it
                }
                else priorityQueue.update_position(neighbor);
                
            }
        }
        if (c_x+1 < h && c_y-1 >= 0)   {                                        // check top left
            DistancePixel* neighbor = &distance_image.GetPixel(c_x+1, c_y-1);
            if (ShouldBeUpdated(currentPixel, neighbor)) {
                
                neighbor->setX(currentPixel->getX()+1);
                neighbor->setY(currentPixel->getY()-1);
                neighbor->setValue(currentPixel->getValue() + sq2);

                if (!priorityQueue.in_heap(neighbor)) {
                    priorityQueue.push(neighbor);
                }
                else priorityQueue.update_position(neighbor);
            }
        }
        if (c_x+1 < h)                 {                                        // check top
            DistancePixel* neighbor = &distance_image.GetPixel(c_x+1, c_y);
            if (ShouldBeUpdated(currentPixel, neighbor)) {
                
                neighbor->setX(currentPixel->getX()+1);
                neighbor->setY(currentPixel->getY());
                neighbor->setValue(currentPixel->getValue() + 1);
                
                if (!priorityQueue.in_heap(neighbor)) {
                    priorityQueue.push(neighbor);
                }
                else priorityQueue.update_position(neighbor);
            }
        }
        if (c_x+1 < h && c_y+1 < w)   {                                         // check top right
            DistancePixel* neighbor = &distance_image.GetPixel(c_x+1, c_y+1);
            if (ShouldBeUpdated(currentPixel, neighbor)) {
                
                neighbor->setX(currentPixel->getX()+1);
                neighbor->setY(currentPixel->getY()+1);
                neighbor->setValue(currentPixel->getValue() + sq2);
                
                if (!priorityQueue.in_heap(neighbor)) {
                    priorityQueue.push(neighbor);
                }
                else priorityQueue.update_position(neighbor);
            }
        }
        if (c_y+1 < w)                 {                                        // check right
            DistancePixel* neighbor = &distance_image.GetPixel(c_x, c_y+1);
            if (ShouldBeUpdated(currentPixel, neighbor)) {
                
                neighbor->setX(currentPixel->getX());
                neighbor->setY(currentPixel->getY()+1);
                neighbor->setValue(currentPixel->getValue() + 1);

                if (neighbor->getValue() > max_distance) {
                    max_distance = neighbor->getValue();
                }
                if (!priorityQueue.in_heap(neighbor)) {
                    priorityQueue.push(neighbor);
                }
                else priorityQueue.update_position(neighbor);
            }
        }
        if (c_y+1 < w && c_x-1 >= 0)    {                                       // check bottom right
            DistancePixel* neighbor = &distance_image.GetPixel(c_x-1, c_y+1);
            if (ShouldBeUpdated(currentPixel, neighbor)) {
                
                neighbor->setX(currentPixel->getX()-1);
                neighbor->setY(currentPixel->getY()+1);
                neighbor->setValue(currentPixel->getValue() + sq2);

                if (!priorityQueue.in_heap(neighbor)) {
                    priorityQueue.push(neighbor);
                }
                else priorityQueue.update_position(neighbor);
            }
        }
        if (c_x-1 >= 0)                 {                                       // check bottom
            DistancePixel* neighbor = &distance_image.GetPixel(c_x-1, c_y);
            if (ShouldBeUpdated(currentPixel, neighbor)) {
                
                neighbor->setX(currentPixel->getX()-1);
                neighbor->setY(currentPixel->getY());
                neighbor->setValue(currentPixel->getValue() + 1);
                
                if (!priorityQueue.in_heap(neighbor)) {
                    priorityQueue.push(neighbor);
                }
                else priorityQueue.update_position(neighbor);
            }
        }
        if (c_x-1 >= 0 && c_y-1 >= 0)   {                                       // check bottom left
            DistancePixel* neighbor = &distance_image.GetPixel(c_x-1, c_y-1);
            if (ShouldBeUpdated(currentPixel, neighbor)) {
                
                neighbor->setX(currentPixel->getX()-1);
                neighbor->setY(currentPixel->getY()-1);
                neighbor->setValue(currentPixel->getValue() + sq2);

                if (!priorityQueue.in_heap(neighbor)) {
                    priorityQueue.push(neighbor);
                }
                else priorityQueue.update_position(neighbor);
            }
        }
        
        priorityQueue.pop();                                                // remove first pixel in priority queue
        
    }
   
    /*
    
    STEPS:
    1. Current pixel == highest priority (top of priority queue)
    2. Check if neighbor exists (boundary checks)
    3. Check if neighbor's D is larger than current D + D2neighbor
    4. Update neighbor in priority queue
    5. Pop current from queue

    */
    
    double max = 0;
    for (int i = 0; i < distance_image.Height(); i++) {
        for (int j = 0; j < distance_image.Width(); j++) {
            if (distance_image.GetPixel(j, i).getValue() > max) {
                max = distance_image.GetPixel(j, i).getValue();
            }
        }
    }
    
    return max;
}


// ===================================================================================================

Color Rainbow(double distance, double max_distance) {
  Color answer;
  if (distance < 0.001) {
    // black
    answer.r = 0; answer.g = 0; answer.b = 0;
  } else if (distance < 0.2*max_distance) {
    // blue -> cyan
    double tmp = distance * 5.0 / max_distance;
    answer.r = 0;
    answer.g = tmp*255;
    answer.b = 255;
  } else if (distance < 0.4*max_distance) {
    // cyan -> green
    double tmp = (distance-0.2*max_distance) * 5.0 / max_distance;
    answer.r = 0;
    answer.g = 255;
    answer.b = (1-tmp*tmp)*255;
  } else if (distance < 0.6*max_distance) {
    // green -> yellow
    double tmp = (distance-0.4*max_distance) * 5.0 / max_distance;
    answer.r = sqrt(tmp)*255;
    answer.g = 255;
    answer.b = 0;
  } else if (distance < 0.8*max_distance) {
    // yellow -> red
    double tmp = (distance-0.6*max_distance) * 5.0 / max_distance;
    answer.r = 255;
    answer.g = (1-tmp*tmp)*255;
    answer.b = 0;
  } else if (distance < max_distance) {
    // red -> white
    double tmp = (distance-0.8*max_distance) * 5.0 / max_distance;
    answer.r = 255;
    answer.g = tmp*255;
    answer.b = tmp*255;
  } else {
    // white
    answer.r = answer.g = answer.b = 255;
  }  
  return answer;
}

Color GreyBands(double distance, double max_value, int num_bands) {
  Color answer;
  if (distance < 0.001) {
    // red
    answer.r = 255; answer.g = 0; answer.b = 0;
  } else {
    // shades of grey
    answer.r = answer.g = answer.b = int(num_bands*256*distance/double(max_value)) % 256;
  }  
  return answer;
}

Color AlienChristmas(double distance, double max_value, int num_bands) {
    Color answer;
    if (distance < 0.001) {
        // minty
        answer.r = 255; answer.g = 10; answer.b = 10;
    } else {
        // shades of blue
        answer.b = int(num_bands*256*distance/double(max_value)) % 256;
        answer.r = rand();
        answer.g = (2*answer.b);
    }
    return answer;
}

// ===================================================================================================
