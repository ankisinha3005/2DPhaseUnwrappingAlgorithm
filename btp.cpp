#include<bits/stdc++.h>

using namespace std;


bool check(double a,double b){                       //A function to check if two double are equal
    if(abs(a-b)<DBL_EPSILON)
       return true; 
    else                                    
       return false;  
}
class phase{
	public:
	double val;
	int row;
	int col;
};
double pi = 3.14159;                                     //pi value

bool comp(const phase a1,const phase a2){               //comparator for sorting the array
	return a1.val<a2.val;
}

double funcwrap1(double x){                            //wrapper
	double y = fmod(x,2*pi)-pi;
	return y;
}
void funcQGPU_Hist(vector<vector<double> > &ph, vector<vector<double> > &phaseQualityMap, vector<vector<int> > &treated, 
                  vector<phase> &sort_values, double thr,vector<vector<double> > &phase2DQG ){
    int x = sort_values[0].row;
    int y = sort_values[0].col;
    if(treated[x][y]==0){
    	phase2DQG[x][y]=ph[x][y];
    	treated[x][y]=1;
	}
	// change n and m 
	int m = ph.size();
	int n = ph[0].size();
    int cnt2=0;
	for(int i=1;i<sort_values.size();i++){
		x = sort_values[i].row;
		y = sort_values[i].col;
		if(treated[x][y]==0){
			double a=2.0,b=2.0,c=2.0,d=2.0,e=2.0,f=2.0,g=2.0,h=2.0;
			if(x>0 && treated[x-1][y]==1 && phaseQualityMap[x-1][y]<=thr)  a = phaseQualityMap[x-1][y];
			if(y>0 && treated[x][y-1]==1 && phaseQualityMap[x][y-1]<=thr)  b = phaseQualityMap[x][y-1];
			if(x<m-1 && treated[x+1][y]==1 && phaseQualityMap[x+1][y]<=thr)c = phaseQualityMap[x+1][y];
			if(y<n-1 && treated[x][y+1]==1 && phaseQualityMap[x][y+1]<=thr)d = phaseQualityMap[x][y+1];
			if(x>0 && y>0 && treated[x-1][y-1]==1 && phaseQualityMap[x-1][y-1]<=thr)  e = phaseQualityMap[x-1][y-1];
			if(x<m-1 && y>0 && treated[x+1][y-1]==1 && phaseQualityMap[x+1][y-1]<=thr)  f = phaseQualityMap[x+1][y-1];
			if(x>0 && y>n-1 && treated[x-1][y+1]==1 && phaseQualityMap[x-1][y+1]<=thr)  g = phaseQualityMap[x-1][y+1];
			if(x<m-1 && y<n-1 && treated[x+1][y+1]==1 && phaseQualityMap[x+1][y+1]<=thr)h = phaseQualityMap[x+1][y+1];
			double val = min(a,min(b,min(c,min(d,min(e,min(f,min(g,h)))))));
			double  temp=DBL_MAX;
               if(val>1)
                continue;
			if(check(val,a))temp = phase2DQG[x-1][y]-ph[x][y];
			if(check(val,b))temp = phase2DQG[x][y-1]-ph[x][y];
			if(check(val,c))temp = phase2DQG[x+1][y]-ph[x][y];
			if(check(val,d))temp = phase2DQG[x][y+1]-ph[x][y];
			if(check(val,e))temp = phase2DQG[x-1][y-1]-ph[x][y];
			if(check(val,f))temp = phase2DQG[x+1][y-1]-ph[x][y];
			if(check(val,g))temp = phase2DQG[x-1][y+1]-ph[x][y];
			if(check(val,h))temp = phase2DQG[x+1][y+1]-ph[x][y];
			
			if(!check(temp,DBL_MAX)){
                temp = 2*pi*round(temp/(2*pi));
                phase2DQG[x][y]=temp+ph[x][y];
                treated[x][y]=1;
                cnt2=1;
            }
			
   		}
       
        if(i==sort_values.size()-1){
            i=0;
            if(cnt2==0)
              break;
            else
              cnt2=0;
        }

	}
                  	
}
vector<vector<double> > funcQualityMap(vector<vector<double> > ph){
	int win = 1;
	int m = ph.size();
	int n = ph[0].size();
	vector<vector<double> > phPad(m+2,vector<double>(n+2,0));

	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++)phPad[i+1][j+1]=ph[i][j];
	}
	//padding the left,right,top and bottom rows with zeroes;
	vector<vector<double> > phaseQualityMap(m+2,vector<double>(n+2,0));
	double maxi = -DBL_MAX;
	
	
	for(int k=win;k<=m-1;k++){
		for(int l=win;l<=n-1;l++){
		double H,V,D1,D2;
		  H = funcwrap1(phPad[k][l-1]- phPad[k][l]) -   funcwrap1(phPad[k][l] - phPad[k][l+1]);
          V = funcwrap1(phPad[k-1][l] - phPad[k][l]) -   funcwrap1(phPad[k][l] - phPad[k+1][l]);
          D1 = funcwrap1(phPad[k-1][l-1] - phPad[k][l]) -   funcwrap1(phPad[k][l] - phPad[k+1][l+1]);
          D2 = funcwrap1(phPad[k+1][l-1] - phPad[k][l]) -   funcwrap1(phPad[k][l] - phPad[k-1][l+1]);
          phaseQualityMap[k][l] = sqrt(H*H + V*V + D1*D1 + D2*D2);
          maxi=max(maxi,phaseQualityMap[k][l]);
	}
 }
	
	vector<vector<double> > phaseQualityMap2(m,vector<double>(n,0));
	
	
	
	
	//reusing the array ph
// 	cout<<"Maxi=="<<maxi<<"\n";
	 for(int i=0;i<m;i++){
		for(int j=0;j<n;j++) {
            phaseQualityMap2[i][j]=(phaseQualityMap[i+1][j+1]/maxi);
        }
	}
	
	for(int i=0;i<m;i++){
		phaseQualityMap2[i][0]=1;
		phaseQualityMap2[i][1]=1;
		phaseQualityMap2[i][m-2]=1;
		phaseQualityMap2[i][m-1]=1;
	}
	for(int i=0;i<n;i++){
		phaseQualityMap2[0][i]=1;
		phaseQualityMap2[1][i]=1;
		phaseQualityMap2[n-1][i]=1;
		phaseQualityMap2[n-2][i]=1;
	}
	return phaseQualityMap2;
}

vector<double> linespace(double min, double max, int n){                    // A function to generate 
                                                                           //equal spaced threshold values
 	vector<double> result;
 	int iterator = 0;
 
	for (int i = 0; i <= n-2; i++){
 		double temp = min + i*(max-min)/(floor((double)n) - 1);
 		result.insert(result.begin() + iterator, temp);
 		iterator += 1;
 	}
    result.insert(result.begin() + iterator, max);
 	return result;
}
void funcQGPU(vector<vector<double> > ph, int nbins){
	vector<vector<double> > phaseQualityMap = funcQualityMap(ph);
 
	int m = ph.size();                                                     //number of rows
	int n = ph[0].size();                                                 //number of columns
	
    vector<phase> sort_values;                                                /*A linear array to store the
                                            	                                phaseQualityMap values and its
                                            	                                locations and sort them according
                                            	                                to the quality value
                                            	                             */ 
	 for(int i=0;i<m;i++){
	     
	 	for(int j=0;j<n;j++){
	 		phase t;
	 		t.val = phaseQualityMap[i][j];
	 		t.row=i;
	 		t.col = j;
	 		sort_values.push_back(t);
		 }
		 
	 }
	 
	 sort(sort_values.begin(),sort_values.end(),comp);  // sorting the array

	 vector<double> histcounts = linespace(0,1,nbins);  //stores the different value of threshold
	 vector<vector<int> > treated (m,vector<int>(n,0)); //which pixels were treated; 
	 vector<vector<double> > phase2DQG (m,vector<double>(n,0));//final answer;
	 
	 //function calling for different values of threshold
	 
	 for(int i=0;i<histcounts.size();i++){
	 	double thr = histcounts[i];
	 	funcQGPU_Hist(ph,phaseQualityMap,treated,sort_values,thr,phase2DQG);
	 }
	 
	 
   // printing the final answer
  for(int i=0;i<10;i++){
        for(int j=0;j<10;j++)
         cout<<phase2DQG[i][j]<<"( "<<treated[i][j]<<")"<<" ";    
         cout<<"\n";
    }


}
int main(){
	int m=10;
	int n=10;
// 	phase = peaks(100) % 100 X 100 matrix
//         S = exp(1i*phase);
//         wrapped_phase = angle(S);

	vector<vector<double> > ph2D(m,vector<double>(n));   //input phase matrix
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			double temp;
			cin>>temp;
			ph2D[i][j]=temp;
        
		}
	}

    funcQGPU(ph2D,10); // Function call in the main function which does the unwrapping
    





	
	
	
	return 0;
	
}

