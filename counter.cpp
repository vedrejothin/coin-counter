/*
RUN:
make
./counter hough in.jpg out.jpg
./counter seg-1 in.jpg out.jpg
./counter seg-2 in.jpg out.jpg
./counter eval out1.jpg out2.jpg
*/

//Link to the header file
#include "CImg.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <limits.h>

//parameters for the 4th question..
#define MIN_VARIANCE 10
#define MAX_VARIANCE 100 // change it to 87 for 6.jpg??
#define MIN_BACK_SIZE 2000 // play with the number of pixels..

#define PI 3.14159265
#define radiusLow 15
#define radiusHigh 35
#define edgeLower 30
#define edgeUpper 80
#define circleThreshold 110

//params for third question

//reducing the below threshold will give us more segments..
#define segmentationThreshold 5
#define k_threshold 5 //TODO: experiment with this value..
#define MIN_CircleSize 100

//below threshold is for checking the minimum number of pixels in the component that must match the color range.
#define MIN_MATCH_COLOR 0.4

#define AspectRatio_MIN 0.85
#define AspectRatio_MAX 1.15

//the fill ratio checks...actual ratio should be 0.785
#define FillRatio_MIN 0.7
#define FillRatio_MAX 0.9


//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

CImg<double> averageGrayscale(CImg<double> input);
CImg<double> compress(CImg<double> input);
CImg<double> edgeDetect(CImg<double> input);
void hysteresis(CImg<double> &input, int i, int j);
int hough(CImg<double> input);



class Node
{
	private:
		double red;
		double blue;
		double green;
	

	public:
		double getRed()
		{
			return red;
		}
	
		double getBlue()
		{
			return blue;
		}
		
		double getGreen()
		{
			return green;
		}
		
		int x,y;
		
		Node()
		{
			//does nothing..
		}
		
		Node(int x, int y, double r, double g, double b) //Note that the order is RGB
		{
			//x and y are the position
			this->x=x;
			this->y=y;
			
			//below three are its values in their respective channels..
			this->red=r;
			this->blue=b;
			this->green=g;	
			currComponent=-1;		
		}
		
		//override the comparator operator
	    bool operator == (const Node& n1) const
    	{
	        if(n1.x==x && n1.y==y)
	        {
	        	return true;
	        }
	        return false;
    	}
		
		int currComponent;				
};


class Edge
{
	private:
	
		void computeWeight()
		{
			weight=0;
			
			//compute red^2
			weight+=pow(incNodes[0]->getRed()-incNodes[1]->getRed(), 2);
			
			//compute blue^2
			weight+=pow(incNodes[0]->getBlue()-incNodes[1]->getBlue(), 2);
			
			//compute green^2
			weight+=pow(incNodes[0]->getGreen()-incNodes[1]->getGreen(), 2);
			
			//get the square root.
			weight=sqrt(weight);				
		}

	public:
	
		//the two nodes incident to the edge..
		Node* incNodes[2];

	
		double weight;
		Edge()
		{
			//empty constructor
			//does nothing..
		}
	
		Edge(Node* n1, Node* n2)
		{			
			incNodes[0]=n1;
			incNodes[1]=n2;
			computeWeight();
		}
	
		double getWeight()
		{
			return weight;
		}
		
		//override the comparator operator
	    bool operator < (const Edge& e1) const
    	{
	        return weight < e1.weight;
    	}
		
};

class Component
{
	public:
		double minX, maxX, minY, maxY;
		int id;
		vector<Node*> liNodes;
		
		double maxEdgeWeight;
		Component(int id, vector<Node*> nodes)
		{
    		//init values..
    		minX=INT_MAX;
    		minY=INT_MAX;
    		maxX=INT_MIN;
    		maxY=INT_MIN;
									
			this->id=id;
			maxEdgeWeight=-1;
			merge(nodes, -1);
		}
		
		void merge(vector<Node*> nodes, double weight)		
		{
			if(weight > maxEdgeWeight)
			{
				maxEdgeWeight=weight;
			}
		
			//add all nodes to this component.
			for(int i=0; i<nodes.size(); i++)
			{
				//change the component status
				nodes[i]->currComponent=id;
				liNodes.push_back(nodes[i]);				
			}		
		}
		
		//override the == operator
	    bool operator == (const Component& c1) const
    	{
	        return id == c1.id;
    	}	
    	
    	int getArea()
    	{
    		//area is just the number of pixels.
    		return liNodes.size();
    	}	
    	
		bool isCircle()
		{
			//checks for the minimum size of the segment..
			if(liNodes.size() < MIN_CircleSize)
			{
				return false;
			}
						
					
			computeProperties();			
			//define aspect ratio as width/height    		    		
			double aRatio = (maxX-minX)/(maxY-minY);
			
			if(aRatio < AspectRatio_MIN || AspectRatio_MAX < aRatio)
			{
				//aspect ratio is too high or too low for it to be a circle.
				//cout<<"Aspect ratio problem"<<aRatio<<endl;
				return false;
			}
			
			//check the fill Ratio..
			double fillRatio = getArea()/((maxX-minX)*(maxY-minY));			
			if( fillRatio < FillRatio_MIN || fillRatio > FillRatio_MAX)
			{

				//Fill ratio is too high or too low for it to be a circle.
				//cout<<"fill ratio problem"<<fillRatio<<endl;
				return false;
			}
			
			//check more properties of the region ??? TODO
			
			return true;
		}
		
		void computeProperties()
		{
			//compute the min and max bounding box
    		for(int i=0; i<liNodes.size(); i++)
    		{
    			//get the boundary pixels..
    		
    			if(liNodes[i]->x < minX)
    			{
    				minX=liNodes[i]->x;
    			}
    			
    			if(liNodes[i]->y < minY)
    			{
    				minY=liNodes[i]->y;
    			}

    			if(liNodes[i]->x > maxX)
    			{
    				maxX=liNodes[i]->x;
    			}

    			if(liNodes[i]->y > maxY)
    			{
    				maxY=liNodes[i]->y;
    			}    			
    		}				
		}
		
		
};


void segmentation(CImg<double> &inp, bool isSimpleSeg, bool removeBack)
{
	vector<Component> liComponent;
	vector< vector<Node> > imageNodes;
	vector<Edge> liEdges;

	for(int i=0; i<inp.width(); i++)
	{
		vector<Node> rowNodes;
		for(int j=0; j<inp.height(); j++)
		{
			//create a new node for every vertex
			Node n(i, j, inp(i, j, 0, 0), inp(i, j, 0, 1), inp(i, j, 0, 2));
			
			rowNodes.push_back(n);		
		}
		imageNodes.push_back(rowNodes);	
	}
	
	//create all edges..
	for(int i=0; i<inp.width(); i++)
	{
		for(int j=0; j<inp.height(); j++)
		{
			//cout<<i<<"\t"<<j<<endl;
			
			//for every vertex add only the edges to its right and to its down..			
			if(i!=inp.width()-1)
			{
				Edge e(&imageNodes[i][j], &imageNodes[i+1][j]);
				if(isSimpleSeg)
				{
				    if(e.weight<segmentationThreshold)
				    {
                        //discard all other edges..<optimization trick>
    				    liEdges.push_back(e);
				    }
				}
				else
				{
				    liEdges.push_back(e);
			    }
			}
			
			if(j!=inp.height()-1)
			{
				Edge e(&imageNodes[i][j], &imageNodes[i][j+1]);
				
				if(isSimpleSeg)
				{
				    if(e.weight<segmentationThreshold)
				    {
                        //discard all other edges..<optimization trick>
    				    liEdges.push_back(e);
				    }
				}
				else
				{
				    liEdges.push_back(e);
			    }
			}
		}
	}	
	
	int compId=0;
	
	//create a new component for every node.
	for(int i=0; i<inp.width(); i++)
	{
		for(int j=0; j<inp.height(); j++)
		{
			vector<Node*> v;
			v.push_back(&imageNodes[i][j]);
			
			//create components with ids..
			Component c(compId, v);
			compId++;			
			
			liComponent.push_back(c);
		}	
	}

	int count = liComponent.size();
	
	//sort all the edges based on edge weight.. (in ascending order)
	std::sort(liEdges.begin(), liEdges.end());
		
	for(int i=0; i<liEdges.size(); i++)
	{	
		//check if the incident nodes belong to the same component or not..
		if(liEdges[i].incNodes[0]->currComponent != liEdges[i].incNodes[1]->currComponent)
		{   
            int id1=liEdges[i].incNodes[0]->currComponent, id2=liEdges[i].incNodes[1]->currComponent;
			Component* it1 = &liComponent[id1];
			Component* it2 = &liComponent[id2];
			
			if(!isSimpleSeg)
			{
				//logic for 3.2
				double val1 = it1->maxEdgeWeight + k_threshold/(double)it1->liNodes.size();
				double val2 = it2->maxEdgeWeight + k_threshold/(double)it2->liNodes.size();
				
				double val=val1;
				
				if(val<val2)
				{
					val=val2;
				}
				
				if(liEdges[i].weight>val)
				{
					continue; //skip those edges which violate this condition..
				}
				
			}
			count--;// merge decreases the count by 1						
			it1->merge(it2->liNodes, liEdges[i].weight);
			
			//delete the memory used by the left out component..<TO REDUCE MEMORY FOOT PRINT>
			vector<Node*> temp;
			it2->liNodes.swap(temp);
		}		
	}
	
	CImg<double> out(inp);
	CImg<double> regions(inp);
	
	int possibleCircles =0;	
	int numCircles=0;
	srand (time(NULL)); //seed random number generator..
	for(int c=0; c<liComponent.size(); c++)
	{
		//generate random numbers between 0 and 255 for the three channels..
		int r=rand()%256;
		int b=rand()%256;
		int g=rand()%256;
		if(removeBack)
    	{	
    	    if( liComponent[c].liNodes.size() > MIN_BACK_SIZE)
		    {
		        //remoe tha back ground here..
    		    for(int l=0; l<liComponent[c].liNodes.size(); l++)
                {
	                int x=liComponent[c].liNodes[l]->x;
	                int y=liComponent[c].liNodes[l]->y;
				
	                inp(x, y, 0, 0) = 255;
	                inp(x, y, 0, 1) = 255;
	                inp(x, y, 0, 2) = 255;			
                }
		    }
		}
		else
		{
		    bool isPotentialCircle=false;
		    if(liComponent[c].isCircle())
		    {
			    /*
			    //draw circle at the specified location..
			    double xCenter = (liComponent[i].minX+liComponent[i].maxX)/2;
			    double yCenter = (liComponent[i].minY+liComponent[i].maxY)/2;
			
			    double radius = (-liComponent[i].minX+liComponent[i].maxX - liComponent[i].minY+liComponent[i].maxY)/2;
			    const unsigned char color[] = {255, 255, 255};
			    out.draw_circle( xCenter, yCenter, radius, color, 1, 0U);		
			    cout<<xCenter<<"\t"<<yCenter<<radius<<endl;
			    */
			    numCircles++;						
			    isPotentialCircle=true;
		    }
				
		    for(int l=0; l<liComponent[c].liNodes.size(); l++)
		    {
			    int x=liComponent[c].liNodes[l]->x;
			    int y=liComponent[c].liNodes[l]->y;
						
			    regions(x, y, 0, 0) = r;
			    regions(x, y, 0, 1) = b;
			    regions(x, y, 0, 2) = g;
						
			    if(isPotentialCircle)
			    {				
				    out(x, y, 0, 0) = r;
				    out(x, y, 0, 1) = b;
				    out(x, y, 0, 2) = g;
			    }	
			    else
			    {
				    out(x, y, 0, 0) = 0;
				    out(x, y, 0, 1) = 0;
				    out(x, y, 0, 2) = 0;
			    }		
		    }
		}
	}
    
    if(!removeBack)
    {
	    cout<<"Total Number of Segments: "<<count<<endl;
	    cout<<"Total Number of Circles: "<<numCircles<<endl;
	    out.save("seg_detected.png");
	    regions.save("regions.png");
    }
}


int main(int argc, char **argv) 
{
	if (argc < 3) {
		cout << "Insufficent number of arguments. Please see documentation" << endl;
		cout << "p3 problemId inputfile" << endl;
		return -1;
	}

	char* inputFile = argv[2];
	CImg<double> input(inputFile);
	
	if (input.spectrum() != 3) {
		cout<<"Please give a RGB image as input"<<endl;
		return -1;
	}

	//input.save("resize.jpg");
	/*double R = 0.0, G = 0.0, B = 0.0;
	int n = 0;
	for(int i = 0; i < input.width(); i++) {
		for(int j = 0; j < input.height(); j++) {
			if(!(input(i, j, 0, 0) == 0 || input(i, j, 0, 1) == 0 || input(i, j, 0, 2) == 0)) {
				n++;
				R += input(i, j, 0, 0);
				G += input(i, j, 0, 1);
				B += input(i, j, 0, 2);
			}
		}
	}
	cout << R/n << " " << G/n << " " << B/n << " " << n << endl; */
						

	if(!strcmp(argv[1], "hough"))
	{	
		input.resize(-20, -20, -100, -100, 5);
		//input.blur(1.5);
		int coins = hough(input);
		cout << "Number of coins: " << coins << endl;
	}
	else if(!strcmp(argv[1], "seg-1"))
	{		
		input.resize(-20, -20, -100, -100, 5);
		cout << "# Finding Circular Regions (Simple Segmentation Algorithm)" << endl;				
		segmentation(input, true, false);		
	}
	else if(!strcmp(argv[1], "seg-2"))
	{		
		input.resize(-20, -20, -100, -100, 5);
		input.blur(1);
	
		cout << "# Finding Circular Regions (Modified Segmentation Algorithm)" << endl;				
		segmentation(input, false, false);
	}
	else if(!strcmp(argv[1], "eval"))
	{
		input.resize(-20, -20, -100, -100, 5);
		//input.blur(1.5);
		int coins = hough(input);
		cout << "Number of coins: " << coins << endl;
	}
	
	return 0;
}


CImg<double> averageGrayscale(CImg<double> input) {
	CImg<double> output(input.width(), input.height(), 1, 1, 0);

	for (int x = 0; x < input.width(); x++)
		for (int y = 0; y < input.height(); y++)
			output(x, y, 0, 0) = 0.3 * input(x, y, 0, 0) + 0.6 * input(x, y, 0, 1) + 0.1 * input(x, y, 0, 2);
	
	return output;
}

CImg<double> edgeDetect(CImg<double> input) {
	int width = input.width();
	int height = input.height();
	
	double H[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};	
	
	CImg<double> Ix(width, height, 1, 1, 0);
	CImg<double> Iy(width, height, 1, 1, 0);
	
	for(int i = 1; i < width - 1; i++){
		for(int j = 1; j < height - 1; j++){
			for(int x = 0; x < 3; x++) {
				for(int y = 0; y < 3; y++) {
					Ix(i, j, 0, 0) += input(i - 1 + x, j - 1 + y, 0, 0) * H[x + y * 3];
					Iy(i, j, 0, 0) += input(i - 1 + x, j - 1 + y, 0, 0) * H[x * 3 + y];
				} 
			}
		}
	}
	
	CImg<double> angle(width, height, 1, 1, 0);
	CImg<double> magnitude(width, height, 1, 1, 0);
	
	for(int i = 1; i < width - 1; i++) {
		for(int j = 1; j < height - 1; j++) {
			magnitude(i, j, 0, 0) = sqrt(Ix(i, j, 0, 0) * Ix(i, j, 0, 0) + Iy(i, j, 0, 0) * Iy(i, j, 0, 0));
			angle(i, j, 0, 0) = atan2(Iy(i, j, 0, 0), Ix(i, j, 0, 0)) * 180 / PI;
		}
	}
	
	CImg<double> output(width, height, 1, 1, 0);
	
	for(int i = 1; i < width - 1; i++) {
		for(int j = 1; j < height - 1; j++) {			
			if((angle(i, j, 0, 0) >= -22.5 && angle(i, j, 0, 0) < 22.5) || (angle(i, j, 0, 0) >= 157.5 && angle(i, j, 0, 0) < -157.5)) {
				if(magnitude(i, j, 0, 0) < magnitude(i, j - 1, 0, 0) || magnitude(i, j, 0, 0) < magnitude(i, j + 1, 0, 0))
					output(i, j, 0, 0) = 0;
				else
					output(i, j, 0, 0) = magnitude(i, j, 0, 0);
			} else if((angle(i, j, 0, 0) >= 22.5 && angle(i, j, 0, 0) < 67.5) || (angle(i, j, 0, 0) >= -157.5 && angle(i, j, 0, 0) < -112.5)) {
				if(magnitude(i, j, 0, 0) < magnitude(i + 1, j + 1, 0, 0) || magnitude(i, j, 0, 0) < magnitude(i - 1, j - 1, 0, 0))
					output(i, j, 0, 0) = 0;
				else
					output(i, j, 0, 0) = magnitude(i, j, 0, 0);
			} else if((angle(i, j, 0, 0) >= 67.5 && angle(i, j, 0, 0) < 112.5) || (angle(i, j, 0, 0) >= -112.5 && angle(i, j, 0, 0) < -67.5)) {
				if(magnitude(i, j, 0, 0) < magnitude(i + 1, j, 0, 0) || magnitude(i, j, 0, 0) < magnitude(i - 1, j, 0, 0))
					output(i, j, 0, 0) = 0;
				else
					output(i, j, 0, 0) = magnitude(i, j, 0, 0);
			} else {
				if(magnitude(i, j, 0, 0) < magnitude(i + 1, j - 1, 0, 0) || magnitude(i, j, 0, 0) < magnitude(i - 1, j + 1, 0, 0))
					output(i, j, 0, 0) = 0;
				else
					output(i, j, 0, 0) = magnitude(i, j, 0, 0);	
			}
		}
	}
	
	for(int i = 0; i < width; i++){
		for(int j = 0; j < height; j++){
			if(output(i, j, 0, 0) != 255) {
				if(output(i, j, 0, 0) > edgeUpper) {
					output(i, j, 0, 0) = 255;
					hysteresis(output, i, j);				
				} else {
					output(i, j, 0, 0) = 0;
				}
			}
		}
	}
	
	return output;
}

void hysteresis(CImg<double> &input, int i, int j) {
	for(int x = 0; x < 3; x++) {
		for(int y = 0; y < 3; y++){
			if(x != 1 && y != 1 && i - 1 + x >= 0 && i - 1 + x < input.width() && j - 1 + y >= 0 && j - 1 + y < input.height() && input(i - 1 + x, j - 1 + y, 0, 0) != 255) {
				if(input(i - 1 + x, j - 1 + y, 0, 0) > edgeLower) {
					input(i - 1 + x, j - 1 + y, 0, 0) = 255;
					hysteresis(input, i - 1 + x, j - 1 + y);				
				} else {
					input(i - 1 + x, j - 1 + y, 0, 0) = 0;
				}
			}
		}
	}
}

int hough(CImg<double> colorInput) {
	CImg<double> input = colorInput.spectrum() == 3 ? averageGrayscale(colorInput) : colorInput;
	
	input.blur(1.5);
	input = edgeDetect(input);
	input.save("edges.png");
	
	int count = 0;	
	
	int width = input.width();
	int height = input.height();
	CImg<double> acc(width, height, radiusHigh - radiusLow + 1, 1, 0);
	const unsigned char red[] = {255, 0, 0}, green[] = {0, 255, 0}, blue[] = {0, 0, 255};
	
	for(int i = 0; i < width; i++) {
		for(int j = 0; j < height; j++) {
			if (input(i, j, 0, 0) == 255) {	
				for(int r = radiusLow; r <= radiusHigh; r++) {
					for(int degree = 0; degree < 360; degree++) {
						double theta = (degree * PI) / 180;
						int cx = (int) round(i - r * cos(theta));
						int cy = (int) round(j - r * sin(theta));
						if(cx >= 0 && cx < width && cy >= 0 && cy < height) {
							acc(cx, cy, r - radiusLow, 0) += 1;
						}					
					}
				}
			}
		}
	}
	
	double max = 0.0, min = 255.0;
	for(int i = 0; i < width; i++) {
		for(int j = 0; j < height; j++) {
			for(int r = radiusLow; r <= radiusHigh; r++) {
				double value = acc(i, j, r - radiusLow, 0);
				if(value > max)
					max = value;
				if(value < min)
					min = value;
			}
		}
	}
	
	CImg<double> output(width, height, 1, 1, 0);
	int * radii =  new int[width * height]; 
	
	for(int i = 0; i < width; i++) {
		for(int j = 0; j < height; j++) {
			int rMax = radiusLow;
			double maxCenter = 0.0;
			for(int r = radiusLow; r <= radiusHigh; r++) {
				if(acc(i, j, r - radiusLow, 0) != 0)
					acc(i, j, r - radiusLow, 0) = (acc(i, j, r - radiusLow, 0) - min) * 255 / (max - min);
				if(acc(i, j, r - radiusLow, 0) > maxCenter) {
					rMax = r;
					maxCenter = acc(i, j, r - radiusLow, 0);					
				}
			}
			output(i, j, 0, 0) = acc(i, j, rMax - radiusLow, 0);
			radii[i + j * width] = rMax;			
		}
	}	
	
	int size = 19;
	int k = (size - 1) / 2;
	CImg<double> colorInput2(colorInput);
	for(int i = k; i < width - k; i++){
		for(int j = k; j < height - k; j++) {
			double localMax = output(i, j, 0, 0);
			int iMax = i, jMax = j;
			for(int x = 0; x < size; x++) {
				for(int y = 0; y < size; y++) {
					if(output(i - k + x, j - k + y, 0, 0) >= localMax) {
						localMax = output(i - k + x, j - k + y, 0, 0);
						iMax = i - k + x;
						jMax = j - k + y;
					} else {
						output(i - k + x, j - k + y, 0, 0) = 0;
					}
				} 
			}
		}
	}
	
	
	
	
	double rSum = 0.0;
	int rCount = 0.0;
	for(int i = k; i < width - k; i++){
		for(int j = k; j < height - k; j++) {
			if(output(i, j, 0, 0) > circleThreshold) {
				rSum += radii[i + j * width];
				rCount++;
				
				colorInput2.draw_circle(i, j, radii[i + j * width], blue, 1, 0U);
			}	
		}
	}  	
	
	colorInput2.save("edges_detected.png");
	double rAvg  = rSum / (double) rCount;
	cout << "AVG: " << rAvg << endl;
	
	
		
	for(int i = k; i < width - k; i++){
		for(int j = k; j < height - k; j++) {
			if(output(i, j, 0, 0) > circleThreshold) {
				//below three values will hold the median of individual channels..				
				double R = 0, G = 0, B = 0;							

				int radii_index = i + j * width;				

				int r = radii[radii_index];
				rSum += r;
				rCount++;
				int n = 0, total = 0;
				int matchColor=0;
				int numWhite=0;
				int matchCoinColor=0;	
								
				for(int x = i - r; x < i + r; x++)
					for(int y = j - r; y < j + r; y++)
					{ 
						if(x>=0 && y>=0 && x<width && y< height && (x - i) * (x - i) + (y - j) * (y - j) <= r * r) 
						{
							n++;
							//dividing by 255..so that the values dont go out of bounds..
							R += colorInput(x, y, 0, 0)/(double)255;
							G += colorInput(x, y, 0, 1)/(double)255;
							B += colorInput(x, y, 0, 2)/(double)255;							
							
							//color checking scheme
							int colorThreshold =5;
							//if(colorInput(x, y, 0, 0) > 60)
							//{
							if(colorInput(x, y, 0, 0) > 50)
							if( colorInput(x, y, 0, 1) >= colorInput(x, y, 0, 2) || (colorInput(x, y, 0, 2)- colorInput(x, y, 0, 1)) < colorThreshold)
							{
								if(colorInput(x, y, 0, 0) >= colorInput(x, y, 0, 1) ||  (colorInput(x, y, 0, 0)- colorInput(x, y, 0, 1)) < colorThreshold)
								{
									matchColor++;
								}															
							}
							
							
							//match with coin's average color..
							//if(colorInput(x, y, 0, 0) > 100 && colorInput(x, y, 0, 0) < 150 && colorInput(x, y, 0, 1) > 30 && colorInput(x, y, 0, 1) < 90 
							//	&& colorInput(x, y, 0, 2) > 20 && colorInput(x, y, 0, 2) < 80) 
							if(colorInput(x, y, 0, 0) > 50)
							{
								matchCoinColor++;
							}			
														
						}
						
					}
				matchCoinColor=n;
				//get the original range of values..		
				R *=255; G *=255; B *=255;		
				R /= n; G /= n; B /= n;

				double variance=0; //variance of individual channels.. 
				for(int x = i - r; x < i + r; x++)
					for(int y = j - r; y < j + r; y++) 
						if( x>=0 && y>=0 && x<width && y< height && (x - i) * (x - i) + (y - j) * (y - j) <= r * r) 
						{
							variance += sqrt(pow(colorInput(x, y, 0, 0)-R, 2) + pow(colorInput(x, y, 0, 1)-G,2)+pow(colorInput(x, y, 0, 2)-B,2));
						}
				variance /=n;
				if(radii[i + j * width] > rAvg - 5 && radii[i + j * width] < rAvg + 5)
				{
					if( matchColor > MIN_MATCH_COLOR*n && matchCoinColor > 0.1*n && variance > MIN_VARIANCE && variance < MAX_VARIANCE) 					
					{
						count++;
						colorInput.draw_circle(i, j, radii[i + j * width], blue, 1, 0U);
						//cout<<"Circle variance is "<<variance<<endl;
					} else {
						//cout<<"Color Match Ratio is "<< (matchCoinColor/(double)n)<<endl;
						//colorInput.draw_circle(i, j, radii[i + j * width], red, 1, 0U);
						//cout << output(i, j, 0, 0) << "R: " << radii[i + j * width] << endl;
						//cout<<"Non-Circle variance is "<<variance<<endl;
						//cout<<"Match Ratio is "<< (matchColor/(double)n)<<endl;
					}
				}
				//cout<<"------------------------------------------------------"<<endl;
			}
		}
	}
	delete [] radii;
	//output.save("acc1.jpg");
	//acc.save("acc.jpg");
	colorInput.save("detected.png");
	return count;
}
