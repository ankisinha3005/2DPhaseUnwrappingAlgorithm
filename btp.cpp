#include<bits/stdc++.h>
//#include "matplotlibcpp.h"
//namespace plt = matplotlibcpp;
using namespace std;

//sudo apt-get install python-matplotlib python-numpy python2.7-dev 
//g++ btp.cpp -std=c++11 -I/usr/include/python2.7 -lpython2.7
bool check(double a,double b){
    if(abs(a-b)<DBL_EPSILON)
      return true;
    return false;  
}
class phase{
	public:
	double val;
	int row;
	int col;
};
double pi = 3.14159;

bool comp(const phase a1,const phase a2){
	return a1.val<a2.val;
}

double funcwrap1(double x){              
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

vector<double> linespace(double min, double max, int n){
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
 
	int m = ph.size();
	int n = ph[0].size();
//   for(int i=0;i<m;i++){
//       for(int j=0;j<n;j++)cout<<phaseQualityMap[i][j]<<" ";
//       cout<<"\n";
//   }

	vector<phase> sort_values; //sorted values of quality map; 
	 for(int i=0;i<m;i++){
	 	for(int j=0;j<n;j++){
	 		phase t;
	 		t.val = phaseQualityMap[i][j];
	 		t.row=i;
	 		t.col = j;
	 		sort_values.push_back(t);
		 }
	 }
	 sort(sort_values.begin(),sort_values.end(),comp);
// 	 for(int i=0;i<sort_values.size();i++)
// 	   cout<<" sort_values ="<<sort_values[i].val<<" i= "<<sort_values[i].row<<" j ="<<sort_values[i].col<<"\n";
	 vector<double> histcounts = linespace(0,1,nbins);
	 vector<vector<int> > treated (m,vector<int>(n,0)); //which pixels were treated; 
	 vector<vector<double> > phase2DQG (m,vector<double>(n,0));//final ans;
	 for(int i=0;i<histcounts.size();i++){
	 	double thr = histcounts[i];
	 //	cout<<histcounts[i]<<" ";
	 	funcQGPU_Hist(ph,phaseQualityMap,treated,sort_values,thr,phase2DQG);
	 }
  
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
	vector<vector<double> > ph2D(m,vector<double>(n));
	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			double temp;
			cin>>temp;
			ph2D[i][j]=temp;
        
		}
	}

    funcQGPU(ph2D,10);
    
   cout<<ph2D[6][7];



//	
//	
//	std::vector<std::vector<double>> x, y, z;
//    for (double i = -5; i <= 5;  i += 0.25) {
//        std::vector<double> x_row, y_row, z_row;
//        for (double j = -5; j <= 5; j += 0.25) {
//            x_row.push_back(i);
//            y_row.push_back(j);
//            z_row.push_back(::std::sin(::std::hypot(i, j)));
//        }
//        x.push_back(x_row);
//        y.push_back(y_row);
//        z.push_back(z_row);
//    }
//
//    plt::plot_surface(x, y, z);
//    plt::show();
	
	
	
	
	
	
	return 0;
	
}

