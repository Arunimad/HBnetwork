#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<math.h>
using namespace std;



// 


//


//
ifstream inputfile;
string word,line;
ifstream file;
//



//
double HBbondlength_cutoff;
double HBbondangle_cutoff;
string A;
string X;
string D;
string E,F;
string dummy;

struct AXD{
string A;
string X;
string D;
}AXD[10];

//
double num_atoms;

struct coordinates {
	string atom;
	double x;
	double y;
	double z;
	int cluster_no =0;

}data_site1[1000000],data_site2[1000000],data_site3[1000000],data_site4[1000000],data_site5[1000000];

int cluster_summed_all_frames[1000][1000];                       // to calculate how many times out of nframes, it stays in a cluster
int cluster_finalised[1000];                                    //  the cluster number it stayed most
//

int temp1=0;
int temp2=0;
int temp3=0;
int temp4=0;
int temp5=0;
int nframes=0;
//


//
int cluster_counter=0;
double dx,dy,dz;
double dist;

double r_pore =10;

double dxx,dyy,dzz;
double distt;

double dxnorm,dynorm,dznorm;
double dxxnorm,dyynorm,dzznorm;

double angle;


//

int main(){


     // Reading from the  INPUT file
	inputfile.open("INPUT");
	getline(inputfile,line);
	inputfile>>word;
	inputfile>>word;
	inputfile>>HBbondlength_cutoff;
	cout<<"\n"<<"HBbondlength_cutoff ="<<HBbondlength_cutoff;
	inputfile>>word;
	inputfile>>word;
	inputfile>>HBbondangle_cutoff;
	cout<<"\n"<<"HBbondangle_cutoff ="<<HBbondangle_cutoff;
	inputfile>>A;
	cout<<"\n"<<A;
	inputfile>>X;
	cout<<"\n"<<X;
	inputfile>>D;
	cout<<"\n"<<D;
	inputfile>>dummy;
	cout<<"\n"<<dummy;
	inputfile.close();


     // 
    
    

     //........................................................................................//


	//Reading HISTORY file
       file.open("HISTORY");


       //skipping first line
       getline(file,line);
       cout<<"\n"<<line;

       //get number of atoms from second line

       file >> word;
       file >> word;
       file >> num_atoms;
       cout<<"\n"<<"No. of atoms = "<<num_atoms;
       getline(file,line);

       //Skipping next 4 lines
       for(int i=0; i<4; i++){
       file.ignore(256,'\n');  
       }

       //store the co-ordinates of AXD in a struct

       while (true){
       
       for(int i=0; i< num_atoms;i++){

	       file>>word;
	       if(word.compare(A)==0){
	        
	        data_site1[temp1].atom = A;       
		getline(file,line);
		file>>data_site1[temp1].x;
		file>>data_site1[temp1].y;
		file>>data_site1[temp1].z;
         	temp1++;
		getline(file,line);
//       	data_site1[temp1].cluster_no = {0};

	       }
	       else if (word.compare(X)==0){
	       
	       data_site2[temp2].atom = X; 
	       getline (file,line);
	       file>>data_site2[temp2].x;
	       file>>data_site2[temp2].y;
	       file>>data_site2[temp2].z;
	       temp2++;
	       getline(file,line);
//	       data_site2[temp2].cluster_no[] = {0};    
	              
	       
	       }
	       else if (word.compare(D)==0){
	       
	       data_site3[temp3].atom = D; 
	       getline (file,line);
	       file>>data_site3[temp3].x;
	       file>>data_site3[temp3].y;
	       file>>data_site3[temp3].z;
	       temp3++;
	       getline(file,line);
//	       data_site3[temp3].cluster_no[] = {0};
	       }
	       else if(word.compare("O_S")==0){
		       
	       data_site4[temp4].atom = E;
	       getline (file,line);
	       file>>data_site4[temp4].x;
	       file>>data_site4[temp4].y;
	       file>>data_site4[temp4].z;
	       temp4++;
	       getline(file,line);
	       }
	       else if(word.compare("H_OH")==0){

	       data_site5[temp5].atom = F;
	       getline (file,line);
	       file>>data_site5[temp5].x;
	       file>>data_site5[temp5].y;
	       file>>data_site5[temp5].z;
	       temp5++;
	       getline(file,line);
	       }
	       else{
		       
	       getline (file,line);
	       getline (file,line);
	       
	       }

       }
       getline(file,line);


       if(!file.eof()){
	       nframes +=1;

	       //skipping 3 lines of cell details
	       for(int i=0; i<3; i++){
		       getline(file,line);    
	             
	       
	       }
       }
       else{
	       nframes +=1;
	       
	       cout<<"\n"<<"nframes = "<<nframes<<"\n";
	       cout<<"\n"<<"Reached end of file!"<<"\n";
	       break;
    
       
       }    
       
       
       
       }
       
       
       cout<<"\n"<<"temp 1  ="<<temp1;
       cout<<"\n"<<"temp 2  ="<<temp2;
       cout<<"\n"<<"temp 3  ="<<temp3;
       cout<<"\n"<<"temp 4  ="<<temp4;
       cout<<"\n"<<"temp 5  ="<<temp5;
       cout<<"\n"<<"nframes ="<<nframes;
       cout<<"\n";
       
       // stored co-ordinates & initialised clusters as zero
       //
       // NOTE : temp values will be from 0 to 1999 , totalling 2000
       




//................................................... MAIN CALCULATION..................................................................//
       //Now calculating length&angles then Checking for hydrogen bond & assigning clusters

       
       
       for(int l=0; l< nframes; l++){                                                          //outer -outer for loop ,looping over all frames

       
       for (int i=l*(temp2/nframes); i<  ( (l+1)*(temp2/nframes) ) ; i++ ){                                          //outer for


       for(int j=l*(temp3/nframes); j <  ( (l+1)*(temp3/nframes) ) ; j++){                                           //inner for


       	      //if (j!=(i+(j*(temp2/nframes)))){ 
	        if(j!= i){
	       
	       // DISTANCE calculation
	       dx = data_site3[j].x-data_site2[i].x;
	       dy = data_site3[j].y-data_site2[i].y;
	       dz = data_site3[j].z-data_site2[i].z;

	       //PBC correction in z-axis, r_pore is assigned 10
	       dz= dz - ((round((dz/ (2.0*r_pore))))*(2.0*r_pore));
	       
	       dist = sqrt ((dx*dx) + (dy*dy) +(dz*dz));
	   //  cout<<"\n"<<"distance "<<i<<"--"<<j<<"="<<dist;


	       //ANGLE calculation

	       dxx = data_site1[i].x-data_site2[i].x;
	       dyy = data_site1[i].y-data_site2[i].y;
	       dzz = data_site1[i].z-data_site2[i].z;

	       dzz= dzz - ((round((dzz/ (2.0*r_pore))))*(2.0*r_pore));

	       distt = sqrt ((dxx*dxx)+(dyy*dyy)+(dzz*dzz));

	       dxnorm = dx/dist;
	       dynorm = dy/dist;
	       dznorm = dz/dist;
	       dxxnorm = dxx/distt;
	       dyynorm = dyy/distt;
	       dzznorm = dzz/distt;

	       angle = acos( (dxnorm*dxxnorm) + (dynorm*dyynorm) + (dznorm*dzznorm) );
	       angle*= 57.324;                                                                //radian to angle
	  //   cout<<"\n"<<"angle"<<i<<"--"<<j<<"="<<angle;




               /* checking the distance & angles --works correct!
	       cout<<"\n"<<data_site3[j].x;
	       cout<<"\n"<<data_site3[j].y;
	       cout<<"\n"<<data_site3[j].z;
	       cout<<"\n"<<data_site2[i].x;
	       cout<<"\n"<<data_site2[i].y;
	       cout<<"\n"<<data_site2[i].z;
	       cout<<"\n"<<data_site1[i].x;
	       cout<<"\n"<<data_site1[i].y;
	       cout<<"\n"<<data_site1[i].z;
	       

	       
	       cout<<"\n"<<"dx="<<dx;
	       cout<<"\n"<<"dy="<<dy;
	       cout<<"\n"<<"dz="<<dz;
	       cout<<"\n"<<"dist="<<dist;

	       cout<<"\n"<<"dxx="<<dxx;
	       cout<<"\n"<<"dyy="<<dyy;
	       cout<<"\n"<<"dzz="<<dzz;
	       cout<<"\n"<<"dist="<<distt;

	       cout<<"\n"<< ((dxnorm*dxxnorm) + (dynorm*dyynorm) + (dznorm*dzznorm));



	       
	       break;

	       */


              // HB cluster assigning

	      
	       
	       if (dist <= HBbondlength_cutoff  && angle <= HBbondangle_cutoff ){

		       if(data_site2[i].cluster_no==0 && data_site3[j].cluster_no==0){
			       cluster_counter++;
			       
			       data_site2[i].cluster_no=cluster_counter;
			       data_site3[i].cluster_no=cluster_counter;
			       data_site3[j].cluster_no=cluster_counter;
			       data_site2[j].cluster_no=cluster_counter;
		       
		       }else if (data_site2[i].cluster_no!=0 && data_site3[j].cluster_no==0){
			       
			       data_site3[j].cluster_no=data_site2[i].cluster_no;
			       data_site2[j].cluster_no=data_site2[i].cluster_no;
		       
		       }else if (data_site2[i].cluster_no==0 && data_site3[j].cluster_no!=0){
			       
			       data_site2[i].cluster_no=data_site3[j].cluster_no;
			       data_site3[i].cluster_no=data_site3[j].cluster_no;
		       
		// THIS PART IS TRICKY - I CALL THIS COALESCENSE OF CLUSTERS - THE POINT WHERE TWO CLUSTERS CLUBS !!	       
		       }else if (data_site2[i].cluster_no!=0 && data_site3[j].cluster_no!=0){
			       int cluster_tobe_changed;
			      
			       if(data_site2[i].cluster_no != data_site3[j].cluster_no){
				       cluster_tobe_changed = data_site3[j].cluster_no;
				       data_site3[j].cluster_no=data_site2[i].cluster_no;
			       }
			      
			       for (int m=l*(temp2/nframes); m <  ( (l+1)*(temp2/nframes) ) ; m++ ){
				        for(int n=l*(temp3/nframes); n <  ( (l+1)*(temp3/nframes) ) ; n++){
						if(data_site3[n].cluster_no == cluster_tobe_changed){
							data_site3[n].cluster_no == data_site2[i].cluster_no;
						}							
					}
					if(data_site2[m].cluster_no == cluster_tobe_changed){
						data_site2[m].cluster_no = data_site2[i].cluster_no;
					}
			       }// end of this double for loop
						
		       }//coalescense-else if bracket
	       
	    
	   
	       
	       }//if
           	      	      
	        
	      }//if(i!=j)

       }//inner for loop

       
       
       }//outer for loop


              //Just to check ------- checking 10 frames out of 2000, just to see how clusters differs
             
              cout<<"\n";
	      cout<<"Frame "<<l+1<<"\n";
	      cout<<"------------------"; 
              for (int i=l*(temp2/nframes); i< (l+1)*(temp2/nframes) ; i++){
	      cout<<"\n"<<"Lies in cluster " <<data_site2[i].cluster_no;
              }

	      cout<<"\n";               

	  

                      // Resetting values to zero , before going to next frame
		             cluster_counter=0;


       }//outer-outer for loop


       ////      CHECKING proximity with MOF framework - H_OH O_S               ////
       
       
       for(int l=0; l< nframes; l++){
       for(int i= l*(temp2/nframes); i< (l+1)*(temp2/nframes) ; i++){
	       for(int p= l*(temp4/nframes); p< (l+1)*(temp4/nframes); p++){

		       dx= data_site4[p].x-data_site2[i].x;
		       dy= data_site4[p].y-data_site2[i].y;
		       dz= data_site4[p].z-data_site2[i].z;

		       dist = sqrt ((dx*dx)+ (dy*dy)+(dz*dz));
		       //cout<<"distance btw mof/imi ="<<dist<<"\n";

		       if(dist <= HBbondlength_cutoff ){
			       //cout<<"coming in"<<"\n";
			       cout<<"distance btw mof/imi ="<<dist<<"\n";
			       if(data_site4[p].cluster_no==0)
				       data_site4[p].cluster_no=data_site2[i].cluster_no;
			       if(data_site4[p].cluster_no!=0 && data_site4[p].cluster_no!=data_site2[i].cluster_no){
			          // data_site2[i].cluster_no=data_site4[p].cluster_no;
			               cout<<"\n"<<"changing"<<" "<<data_site2[i].cluster_no<<" to "<<data_site4[p].cluster_no;
				       data_site2[i].cluster_no=data_site4[p].cluster_no; }

		       }

	       }//inner for
       
       }//outer for
       }//outer-outer



       
       for(int l=0; l< nframes; l++){
       for(int j= l*(temp3/nframes); j< (l+1)*(temp3/nframes) ; j++){
	       for(int q= l*(temp5/nframes); q< (l+1)*(temp5/nframes); q++){

		      // cout<<"\n"<<"j= "<<j<<" ;; "<<"q= "<<q;
		       dx= data_site5[q].x-data_site3[j].x;
		       dy= data_site5[q].y-data_site3[j].y;
		       dz= data_site5[q].z-data_site3[j].z;

		       dist = sqrt ((dx*dx)+ (dy*dy)+(dz*dz));

		       if(dist <= HBbondlength_cutoff ){

			       if(data_site5[q].cluster_no==0){
				//       cout<<"\n"<<"changing"<<" "<<data_site5[q].cluster_no<<" to "<<data_site2[j].cluster_no;
				       data_site5[q].cluster_no=data_site2[j].cluster_no;
			       }
			       else{
				       cout<<"coming in";
				       cout<<"\n"<<"changing"<<" "<<data_site2[j].cluster_no<<" to "<<data_site5[q].cluster_no;
				       data_site2[j].cluster_no=data_site5[q].cluster_no;
			       }

		       }

	       }//inner for
       
       }//outer for
       }//outer-outer
















       ////        AVERAGING cluster numbers by considering all frames           ////
       
       // Keeping l=1, instead of l=0, to avoid initial frame, that might be oddly config sometimes


	       for (int i=(temp2/nframes); i<  temp2 ; i++ ){     //inner for
		       cluster_summed_all_frames[i%(temp2/nframes)][data_site2[i].cluster_no]++;
	       }//closing inner for loop

       







       //DISPLAY clusters added over all frames
       cout<<"\n"<<"FREQUENCY OVER NFRAMES"<<"\n";
       cout<<"\n"<<"----------------------";
	       for(int i=0;i < (temp2/nframes); i++){
	       for(int j=0; j<10; j++){
		       cout<<"\n"<<i<<"--"<<j<<" "<<cluster_summed_all_frames[i][j];
	       
	       }
       }

	       
       //finalising cluster no. by calculating where it stayed most out of all frames
               int maximum_stay =0;
	       int maximum_stay_clusterno =0;
	       for(int i=0; i<(temp2/nframes); i++){
		       for(int j=0; j< 10; j++){
			       if(cluster_summed_all_frames[i][j] >= maximum_stay){
				       maximum_stay=cluster_summed_all_frames[i][j];
				       maximum_stay_clusterno=j;				       			       
			       }
			       
		       }
		       cluster_finalised[i]=maximum_stay_clusterno;
		       maximum_stay =0;
		       maximum_stay_clusterno=0;
	       }
	       cout<<"\n"<<"\n";
	  


       //Displaying the cluster no. where the molecule stayed most
       
	       for(int i=0; i< (temp2/nframes); i++){
		       cout<<"molecule"<<"-"<<i<<" "<<"highly stayed in cluster no.= "<<cluster_finalised[i];
		       cout<<"\n";
	       }



 



       
       
}// int main()       

