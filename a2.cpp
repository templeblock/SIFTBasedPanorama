// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "Sift.h"
#include <time.h>
typedef pair<string,int> Pair;

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

/*
Transform image using inverse warping
*/
CImg<unsigned char> transform_image(CImg<unsigned char> &input_image, CImg<double> &M,
											int minx=0, int miny=0, int maxx=0, int maxy=0)
{	
	if (minx == 0 && miny == 0 && maxx == 0 && maxy == 0)
	{
		maxx = input_image.width()-1;
		maxy = input_image.height()-1;
	}

	CImg<unsigned char> output_image(maxx-minx+1, maxy-miny+1, 1, 3, 0);
	double w;
	int xsource, ysource;

	M.invert(true);

	for (int x = minx; x <= maxx; ++x)
	{
		for(int y = miny; y <= maxy; ++y)
		{
			w = x*M(0,2)+y*M(1,2)+M(2,2);
			xsource = (int) ((x*M(0,0)+y*M(1,0)+M(2,0)) / w + 0.5);
			ysource = (int) ((x*M(0,1)+y*M(1,1)+M(2,1)) / w + 0.5);

			for (int p = 0; p < 3; ++p)			
				if (!(xsource < 0 || ysource < 0 || xsource >= input_image.width() || ysource >= input_image.height()))
					output_image(x-minx, y-miny, 0, p) = input_image(xsource, ysource, 0, p);			
		}
	}
	return output_image;
}

/*
Find best transformation from point2 to point1 using RANSAC method
threshold - max pixel distance allowed in RANSAC transformation
*/
CImg<double> ransac(vector< pair< pair<int,int>, pair<int,int> > > &matches, int steps, double threshold)
{
	int n = matches.size(), index, max_inliers = -1;
	int x, y, xp, yp;	
	double x1, y1, w;

	CImg<double> A(8, 8, 1, 1, 0), B(1, 8, 1, 1, 0), T(3, 3, 1, 1, 0), T_best;

	time_t t;
	srand((unsigned) time(&t));

	// do this #steps times
	while (steps--)
	{

		// Building the matrices for four points
		for (int i = 0; i < 4; ++i)
		{
			// Getting four points at random
			index = rand() % n;
			xp = matches[index].first.first;
			yp = matches[index].first.second;
			x = matches[index].second.first;
			y = matches[index].second.second;

			// Intialize A and B
			A(0,2*i) = x;	A(1,2*i) = y;	A(2,2*i) = 1;	A(6,2*i) = -x*xp;	A(7,2*i) = -y*xp;
			A(3,2*i+1) = x;	A(4,2*i+1) = y;	A(5,2*i+1) = 1;	A(6,2*i+1) = -x*yp;	A(7,2*i+1) = -y*yp;

			B(0,2*i) = xp;
			B(0,2*i+1) = yp;		
		}

		// Solve 8 linear equations with 8 unknown to get Transform value into B
		B.solve(A);

		// Get the B (1x8) matrix values into the tranform matrix T (3x3)
		for (int i = 0; i < 8; ++i)		
			T(i%3, i/3) = B(0,i);
		T(2,2) = 1;

		// Count number of inliers for this transformation
		int inliers = 0;		
		
		for (int i = 0; i < n; ++i)
		{
			xp = matches[i].first.first;
			yp = matches[i].first.second;
			x = matches[i].second.first;
			y = matches[i].second.second;

			w = x*T(0,2)+y*T(1,2)+T(2,2);
			x1 = (x*T(0,0)+y*T(1,0)+T(2,0)) / w;
			y1 = (x*T(0,1)+y*T(1,1)+T(2,1)) / w;

			if (abs(x1-xp) < threshold && abs(y1-yp) < threshold) // This point matches
				++inliers;			
		}

		// if this projection is current best projection
		if (inliers > max_inliers)
		{
			max_inliers = inliers;
			T_best = T;
		}
	}

	cout << "RANSAC: inlier found " << max_inliers << " out of " << n << "." << endl;

	return T_best;
}

/*
Find sift matches between two images
threshold1 - max euclidian distance allowed
threshold2 - minimum ratio*10.0 between closest match and second closest match allowed
*/
vector< pair< pair<int,int>, pair<int,int> > > get_sift_matches(CImg<unsigned char> &image1, CImg<unsigned char> &image2,
																int threshold1, double threshold2)
{
	vector< pair< pair<int,int>, pair<int,int> > > matches;	

	CImg<double> gray1 = image1.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> descriptors1 = Sift::compute_sift(gray1);

	CImg<double> gray2 = image2.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray2);

	vector<vector<int> > flag_matrix(descriptors1.size(), vector<int>(descriptors2.size(), -1));

	int best_match_index;
	int best_match_value, second_best_match_value;
	
	for(int i = 0; i < descriptors1.size(); i++)
	{
		best_match_value = threshold1*threshold1;
		second_best_match_value = threshold1*threshold1;
		best_match_index = -1;			

		for(int j = 0; j < descriptors2.size(); j++)
		{
			int distance = 0;
			for(int l = 0; l < 128; l++)
			{
				int d;
				d = descriptors1[i].descriptor[l] - descriptors2[j].descriptor[l];				
				distance += d * d;

				if (distance > second_best_match_value)
					break;
			}			

			if (distance < best_match_value)
			{
				if (best_match_index != -1)				
					second_best_match_value = best_match_value;									

				best_match_value = distance;
				best_match_index = j;
			}
			else if (distance < second_best_match_value)
			{
				second_best_match_value = distance;				
			}
		}

		if (best_match_index != -1)
		{
			if (10 * 10 * best_match_value < threshold2 * threshold2 * second_best_match_value)
			{
				flag_matrix[i][best_match_index] = best_match_value;						
			}
		}
	}


	for(int i = 0; i < descriptors2.size(); i++)
	{
		best_match_value = threshold1*threshold1;		
		best_match_index = -1;			

		for(int j = 0; j < descriptors1.size(); j++)
		{
			if (flag_matrix[j][i] == -1)
				continue;

			int distance = flag_matrix[j][i];

			if (distance < best_match_value)
			{				
				best_match_value = distance;
				best_match_index = j;
			}
		}

		if (best_match_index != -1)
		{
			int x1, y1, x2, y2;
			x1 = descriptors1[best_match_index].col;
			y1 = descriptors1[best_match_index].row;
			x2 = descriptors2[i].col;
			y2 = descriptors2[i].row;

			matches.push_back(make_pair(make_pair(x1,y1),make_pair(x2,y2)));
		}
	}

	return matches;
		
}

/*
Create a warped and panorama by combining image1 and image2
thresh1 - max euclidian distance allowed is SIFT matching
thresh2 - minimum ratio*10.0 between closest match and second closest match allowed in SIFT matching
thresh3 - max pixel distance allowed in RANSAC transformation
*/
pair<CImg<unsigned char>,CImg<unsigned char> >create_warped_panorama(CImg<unsigned char> &image1, CImg<unsigned char> &image2,
											int thresh1=200, double thresh2=9.0, double thresh3=10.0)
{
	double x1, y1, w;
	int x, y;
	vector< pair< pair<int,int>, pair<int,int> > > matches;

	matches = get_sift_matches(image1, image2, thresh1, thresh2);	

	cout << "SIFT: " << matches.size() << " matches found." << endl;
	
	CImg<unsigned char> image3(image1.width()+image2.width(), max(image1.height(), image2.height()), 1, 3, 0);
	for (int x = 0; x < image1.width(); ++x)
		for (int y = 0; y < image1.height(); ++y)
			for(int p=0; p<3; p++)
				image3(x, y, 0, p) = image1(x, y, 0, p);

	for (int x = 0; x < image2.width(); ++x)
		for (int y = 0; y < image2.height(); ++y)
			for(int p=0; p<3; p++)
				image3(x+image1.width(), y, 0, p) = image2(x, y, 0, p);
	
	const unsigned char color[] = { 255, 255, 0};
	for (int i = 0; i < matches.size(); ++i)
	{
 		image3.draw_line(matches[i].first.first, matches[i].first.second, matches[i].second.first+image1.width(), matches[i].second.second, color);
	}
	image3.save("image3.png");	
	
	CImg<double> T = ransac(matches, 1000, thresh3);

	// Find the boundary we will get after stitching two images
	int minx = 0, miny = 0, maxx = image1.width(), maxy = image1.height();

	// Check for 4 corner points of second image after translation
	for (int i = 0; i < 4; ++i)
	{
		x1 = i%2?image2.width():0;
		y1 = i/2?image2.height():0;
		
		w = x1*T(0,2)+y1*T(1,2)+T(2,2);
		x = (x1*T(0,0)+y1*T(1,0)+T(2,0)) / w + 0.5;
		y = (x1*T(0,1)+y1*T(1,1)+T(2,1)) / w + 0.5;
		if (x < minx) minx = x;
		if (x > maxx) maxx = x;
		if (y < miny) miny = y;
		if (y > maxy) maxy = y;
	}		

	CImg<unsigned char> warped_image = transform_image(image2, T, minx, miny, maxx, maxy);
	CImg<unsigned char> panorama_image = warped_image;
	for (int x = 0; x < image1.width(); ++x)
		for (int y = 0; y < image1.height(); ++y)
		{
			if (panorama_image(x-minx, y-miny, 0, 0) + panorama_image(x-minx, y-miny, 0, 1) + panorama_image(x-minx, y-miny, 0, 2) == 0)
			{
				for(int p=0; p<3; p++)				
					panorama_image(x-minx, y-miny, 0, p) = image1(x, y, 0, p);
			}
			else
			{
				for(int p=0; p<3; p++)						
					panorama_image(x-minx, y-miny, 0, p) = (panorama_image(x-minx, y-miny, 0, p) + image1(x, y, 0, p) ) / 2;
			}	
		}

	return make_pair(warped_image, panorama_image);
}

int SIFT_match(CImg<double> img1, CImg<double> img2, float threshold2, int question)
{
int threshold1 = 200;
int counter = 0;
int x, y, xp, yp;	

vector< pair< pair<int,int>, pair<int,int> > > matches;	

	CImg<double> gray1 = img1.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> descriptors_1 = Sift::compute_sift(gray1);

	CImg<double> gray2 = img2.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> descriptors_2 = Sift::compute_sift(gray2);

	vector<vector<int> > flag_matrix(descriptors_1.size(), vector<int>(descriptors_2.size(), 0));
	int best_match_index;
	int best_match_value, second_best_match_value;
	
	for(int i = 0; i < descriptors_1.size(); i++)
	{
		best_match_value = threshold1*threshold1;
		second_best_match_value = threshold1*threshold1;
		best_match_index = -1;			

		for(int j = 0; j < descriptors_2.size(); j++)
		{
			int distance = 0;
			for(int l = 0; l < 128; l++)
			{
				int d;
				d = descriptors_1[i].descriptor[l] - descriptors_2[j].descriptor[l];				
				distance += d * d;

				if (distance > second_best_match_value)
					break;
			}			

			if (distance < best_match_value)
			{
				if (best_match_index != -1)				
					second_best_match_value = best_match_value;									

				best_match_value = distance;
				best_match_index = j;
				
			}
			else if (distance < second_best_match_value)
			{
				second_best_match_value = distance;				
			}
		}

		if (best_match_index != -1)
		{
			if (10 * 10 * best_match_value < threshold2 * threshold2 * second_best_match_value)
			{
				flag_matrix[i][best_match_index] = 1;						
			}
		}
	}


	for(int i = 0; i < descriptors_2.size(); i++)
	{
		best_match_value = threshold1*threshold1;		
		best_match_index = -1;			

		for(int j = 0; j < descriptors_1.size(); j++)
		{
			if (flag_matrix[j][i] == 0)
				continue;

			int distance = 0;
			for(int l = 0; l < 128; l++)
			{
				int d;
				d = descriptors_2[i].descriptor[l] - descriptors_1[j].descriptor[l];
				distance += d * d;
			}			

			if (distance < best_match_value)
			{				
				best_match_value = distance;
				best_match_index = j;
			}
		}

		if (best_match_index != -1)
		{
			int x1, y1, x2, y2;
			x1 = descriptors_1[best_match_index].col;
			y1 = descriptors_1[best_match_index].row;
			x2 = descriptors_2[i].col;
			y2 = descriptors_2[i].row;

			matches.push_back(make_pair(make_pair(x1,y1),make_pair(x2,y2)));
		}
	}
	
	CImg<unsigned char> image3(img1.width()+img2.width(), max(img1.height(), img2.height()), 1, 3, 0);
	for (int x = 0; x < img1.width(); ++x)
		for (int y = 0; y < img1.height(); ++y)
			for(int p=0; p<3; p++)
				image3(x, y, 0, p) = img1(x, y, 0, p);

	for (int x = 0; x < img2.width(); ++x)
		for (int y = 0; y < img2.height(); ++y)
			for(int p=0; p<3; p++)
				image3(x+img1.width(), y, 0, p) = img2(x, y, 0, p);
	
	const unsigned char color[] = { 255, 255, 0};
	for (int i = 0; i < matches.size(); ++i)
	{
		counter++;
 		image3.draw_line(matches[i].first.first, matches[i].first.second, matches[i].second.first+img1.width(), matches[i].second.second, color);
	}
	if(question == 1){
	image3.save("part1_normal.png");
	cout << "The merged image has been saved as part1_normal.png in the root folder" << endl;
	}
	return counter;
}

bool sort_vectors(pair<string, int> p1, pair<string, int> p2)
{	try{
       if(p1.second <= p2.second) return false; 
       else return true;
       }
    catch(const string &err) {
    	cerr << "Error: " << err << endl;
    	cout << "Oops! Don't know what went wrong";
    }
}

//This function is to sort the images with maximum matching points and to find the precision

void maximum_key_points(int argc, char** argv, CImg<double> input_image, string inputFile, int user_input){
    vector<Pair> pair(argc-3);
    int limit_for_comparison = 10;
    int correct_prediction = 0;
    int false_prediction = 0;
	int num_key_points;
	for(int counter = 3 ; counter < argc; counter++)
	{
			CImg<double> img2(string(argv[counter]).c_str());
			num_key_points = SIFT_match(input_image,img2, 9.0, user_input);
			pair.push_back(make_pair(string(argv[counter]).c_str(),num_key_points));
			cout << num_key_points << " feature points matched with source file "<< inputFile << " and target file " << argv[counter] << endl;
	}
	
	sort(pair.begin(), pair.end(), sort_vectors);

	
	if(user_input == 2){
	cout << "-------------------------------------------------------------------------------------------------------" << endl;
		for(int i = 0; i < pair.size() ; i++)
		{
			if(!pair[i].first.empty())
	     	cout << pair[i].first << " has " << pair[i].second << " feature point matches with " << inputFile << endl;
		}
	cout << "-------------------------------------------------------------------------------------------------------" << endl;
	}
	
	else if(user_input == 3){
		cout << "Entered userinput3"<< endl;
		cout << "-------------------------------------------------------------------------------------------------------" << endl;
		cout << "				Top 10 matches for " << inputFile << "								" << endl;
		cout << "-------------------------------------------------------------------------------------------------------" << endl;
		
		for(int i = 0; i < limit_for_comparison ; i++){
			cout << pair[i].first << " has " << pair[i].second << " feature point matches with " << inputFile << endl;
		}
		cout << "-------------------------------------------------------------------------------------------------------" << endl;
		cout << "				Precision calculation for " << inputFile << "						" << endl;
		cout << "-------------------------------------------------------------------------------------------------------" << endl;
		
		for(int i = 0; i < limit_for_comparison ; i++){
			size_t index = inputFile.find_last_of("_"); 
			string just_name = inputFile.substr(0, index);

			if (pair[i].first.find(just_name) != std::string::npos) {
    			correct_prediction+=1;
			}
			else{
				false_prediction+=1;
			}
			if(!pair[i].first.empty())
			cout << pair[i].first << " has " << pair[i].second << " feature point matches with " << inputFile << endl;
		}
		
		float prec = ((correct_prediction * 1.0) / (correct_prediction + false_prediction) * 100.0) ;
		cout << "-------------------------------------------------------------------------------------------------------" << endl;
		cout << "					 The precision is " << prec << "%						"<< endl ;
		cout << "-------------------------------------------------------------------------------------------------------" << endl;
	
	}
}

double gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
    phase = 1 - phase;
    return X;
}

vector<double> SIFT_summary(vector< vector<double> > x,vector<float> v, int k, double w)
{
    vector<double> output;
    for(int i=0;i<k;i++)
        {
        double summary = 0.0;
        for(int j=0;j<128;j++)
            summary +=x[i][j]*v[j];
        output.push_back(summary/w);
        }
    return output;
}

int SIFT_summary_match(CImg<double> img1, CImg<double> img2, int k, int w)
{
	int output = 0;
    int candis = 10;
    vector< vector<double> > x;
    for(int i=0;i<k;i++)
        {
        vector<double> y;
        for(int j=0;j<128;j++)
            y.push_back(gaussrand());
        x.push_back(y);
        }   
	CImg<double> gray_img1 = img1.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> descriptors_1 = Sift::compute_sift(gray_img1);
	CImg<double> gray_img2 = img2.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> descriptors_2 = Sift::compute_sift(gray_img2);
    vector<vector<double> > img1_summary, img2_summary;
    for(int i=0; i<descriptors_1.size();i++)
    {
        vector<double> summary = SIFT_summary(x, descriptors_1[i].descriptor, k, w);
        img1_summary.push_back(summary);
    }  
    for(int i=0; i<descriptors_2.size();i++)
    {   vector<double> summary = SIFT_summary(x, descriptors_2[i].descriptor, k, w);
        img2_summary.push_back(summary);
    }
    vector<vector<double> > distance;
    for(int i=0; i<descriptors_1.size(); i++)
    {
        vector<double> one_row_dis;
        for(int j=0; j<descriptors_2.size(); j++)
        {
            int dis = 0;
            for(int kk=0; kk<k; kk++)
                { 
                    double a = img1_summary[i][kk]-img2_summary[j][kk];
                    dis += sqrt(a*a);
                }
            one_row_dis.push_back(dis);
        }
        distance.push_back(one_row_dis);
    }

    vector<int> closests;
    for(int i=0; i<distance.size();i++)
    {
        double closest = 10000000;
        double second =  10000000;
        int n1, n2;
        for(int j=0; j<distance[0].size();j++)
        {
            if(distance[i][j]<second)
            {
                if(distance[i][j]<closest)
                {
                second = closest;
                closest = distance[i][j];
                n2 = n1;
                n1 = j;
                }
                else
                {
                second = distance[i][j];
                n2 = j;
                }
            }
        }
        closests.push_back(n1);
        closests.push_back(n2);
     } 
    char color[]  = {0,255,0};
    int height = img1.height();
    if(img1.height()<img2.height())
        height = img2.height();

	CImg<double> img3(img1.width()+img2.width(), height,1,3,0);
    for(int i=0; i<img1.width(); i++)
    for(int j=0; j<img1.height(); j++)
    for(int k=0; k<3; k++)
    {
    img3(i,j,0,k) = img1(i,j,0,k);
    }
    for(int i=0; i<img2.width(); i++)
    for(int j=0; j<img2.height(); j++)
    for(int k=0; k<3; k++)
    img3(i+img1.width(),j,0,k) = img2(i,j,0,k);
    for(int i=0; i<descriptors_1.size(); i++)
        {
            double e_dis1=0, e_dis2=0;
            for(int j=0;j<128;j++)
                {

                    double a = descriptors_1[i].descriptor[j] - descriptors_2[closests[i*2]].descriptor[j];
                    e_dis1 += sqrt(a*a);
                    
                    double b = descriptors_1[i].descriptor[j] - descriptors_2[closests[i*2+1]].descriptor[j];
                    e_dis2 += sqrt(b*b);
                } 
            double ratio = e_dis1/e_dis2;              
            int index = closests[i*2];
            if(e_dis2<e_dis1)
                {
                ratio = e_dis2/e_dis1;
                index = closests[i*2+1];
                }
            if(ratio<0.4)       
            {
        img3.draw_line(descriptors_1[i].col,descriptors_1[i].row,descriptors_2[index].col+img1.width(),descriptors_2[index].row,color,1);
        output++;
        }
        }      
	img3.save("part1_fast.png");
	cout << "The merged image has been saved as part1_fast.png in the root folder" << endl;
	return output;
}



int main(int argc, char **argv)
{
  try {
  
    if(argc < 2)
    {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    a2-p1 part_id ..." << endl;
	return -1;
    }
     
     //Initializing the variables
     
     string part = argv[1];
     string inputFile = argv[2];
     string inputFile_2 = argv[3];
     int user_input;
     int part1_question;
     CImg<double> input_image(inputFile.c_str());
	 CImg<double> input_image_2(inputFile_2.c_str());  

    if(part == "part1")
    {
		//Part 1 contains the solution for all the part1 questions.
		cout << "You have entered part 1. Please Enter for which question you need to find the solution" << endl;
		cout <<"The options are as below, " << endl << "1. Press 1 for Question 1" << endl << "2. Press 2 for Question 2" << endl << "3. Press 3 for Question 3" << endl;
		cin >> part1_question;
	
		if(part1_question == 1){
    		int match_numbers = SIFT_match(input_image, input_image_2, 9.0, part1_question);
    	}
    
    	if(part1_question == 2 || part1_question == 3){
			maximum_key_points(argc, argv, input_image, inputFile, part1_question);
		}
	
    }
    
    else if(part == "part1fast"){
    
    //part1fast will contain the solution for 4th question.
    int thresh_fast = 10;
    int thresh_fast1 = 260;
    int match_numbers = SIFT_summary_match(input_image, input_image_2, thresh_fast, thresh_fast1);
    cout << "The total number of matching points are : " << match_numbers;
    }
    
    
    else if(part == "part2")
      {
		CImg<unsigned char> image1(argv[2]);

			for (int i = 3; i < argc; ++i)
			{
				CImg<unsigned char> image2(argv[i]);
				cout << "Processing image: " << argv[i] << endl;
				pair<CImg<unsigned char>, CImg<unsigned char> > warped_panorama = create_warped_panorama(image1, image2);
				string output_name = string(argv[i]);
				output_name.insert(output_name.find("."), "_warped");					
				(warped_panorama.first).save(output_name.c_str());
				cout << "Generated warped image: "	<< output_name << endl;
				output_name = string(argv[i]);
				output_name.insert(output_name.find("."), "_panorama");				
				(warped_panorama.second).save(output_name.c_str());
				cout << "Generated panorama image: " << output_name << endl;
				cout << endl;
			}
      }
    else
      throw std::string("unknown part!");

    // feel free to add more conditions for other parts (e.g. more specific)
    //  parts, for debugging, etc.
  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}








