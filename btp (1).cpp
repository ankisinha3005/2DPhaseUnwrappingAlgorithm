#include<bits/stdc++.h>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std;
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
void funcQGPU_Hist(vector<vector<double> >ph, vector<vector<double> > &phaseQualityMap, vector<vector<int> > &treated, 
                  vector<phase>sort_values, int thr,vector<vector<double> > &phase2DQG ){
    int x = sort_values[0].row;
    int y = sort_values[0].col;
    if(treated[x][y]==0){
    	phase2DQG[x][y]=ph[x][y];
    	treated[x][y]=1;
	}
	int n = ph.size();
	int m = ph[0].size();
	for(int i=0;i<sort_values.size();i++){
		x = sort_values[i].row;
		y = sort_values[i].col;
		if(treated[x][y]==1){
			int a,b,c,d,e,f,g,h;
			a=2;
			b=2;
			c=2;
			d=2;
			e=2;
			f=2;
			g=2;
			if(x>0 && treated[x-1][y]==1 && phaseQualityMap[x-1][y]<thr)  a = phaseQualityMap[x-1][y];
			if(y>0 && treated[x][y-1]==1 && phaseQualityMap[x][y-1]<thr)  b = phaseQualityMap[x][y-1];
			if(x<n-1 && treated[x+1][y]==1 && phaseQualityMap[x+1][y]<thr)c = phaseQualityMap[x+1][y];
			if(y<m-1 && treated[x][y+1]==1 && phaseQualityMap[x][y+1]<thr)d = phaseQualityMap[x][y+1];
			if(y>0 && x>0 && treated[x-1][y-1]==1 && phaseQualityMap[x-1][y-1]<thr)  e = phaseQualityMap[x-1][y-1];
			if(y>0 && x<n-1 && treated[x][y-1]==1 && phaseQualityMap[x+1][y-1]<thr)  f = phaseQualityMap[x+1][y-1];
			if(y<m-1 && x>0 && treated[x][y-1]==1 && phaseQualityMap[x-1][y+1]<thr)  g = phaseQualityMap[x-1][y+1];
			if(y<m-1 && x<n-1 && treated[x][y-1]==1 && phaseQualityMap[x+1][y+1]<thr)h = phaseQualityMap[x+1][y+1];
			int val = min(a,min(b,min(c,min(d,min(e,min(f,min(g,h)))))));
			int temp=INT_MAX;
			if(val==a && a!=2)temp = phase2DQG[x-1][y]-ph[x][y];
			if(val==b && b!=2)temp = phase2DQG[x][y-1]-ph[x][y];
			if(val==c && c!=2)temp = phase2DQG[x+1][y]-ph[x][y];
			if(val==d && d!=2)temp = phase2DQG[x][y+1]-ph[x][y];
			if(val==e && e!=2)temp = phase2DQG[x-1][y-1]-ph[x][y];
			if(val==f && f!=2)temp = phase2DQG[x+1][y-1]-ph[x][y];
			if(val==g && g!=2)temp = phase2DQG[x-1][y+1]-ph[x][y];
			if(val==h && h!=2)temp = phase2DQG[x+1][y+1]-ph[x][y];
			if(temp!=INT_MAX){
		    temp = 2*pi*round(temp/(2*pi));
			phase2DQG[x][y]=temp+ph[x][y];
			treated[x][y]=1;}
			
		}
	}
                  	
}
vector<vector<double> > funcQualityMap(vector<vector<double> > ph){
	int win = 1;
	int m = ph.size();
	int n = ph[0].size();
	vector<vector<double> > phPad(m+2,vector<double>(n+2,0));
	for(int i=1;i<=m-1;i++){
		for(int j=1;j<=n-1;j++) phPad[i][j]=ph[i-1][j-1];
	} //padding the left,right,top and bottom rows with zeroes;
	vector<vector<double> > phaseQualityMap(m+2,vector<double>(n+2,0));
	double max = INT_MIN;
	//below loop
	// i think there is  a problem in the loop ,it should run for  i=1 : i<m-1
	// similarly for j =1 : j<n-1
	
	for(int k=win;k<=m;k++){
		for(int l=win;l<=n;l++){
		double H,V,D1,D2;
		  H = funcwrap1(phPad[k][l-1]- phPad[k][l]) -   funcwrap1(phPad[k][l] - phPad[k][l+1]);
          V = funcwrap1(phPad[k-1][l] - phPad[k][l]) -   funcwrap1(phPad[k][l] - phPad[k+1][l]);
          D1 = funcwrap1(phPad[k-1][l-1] - phPad[k][l]) -   funcwrap1(phPad[k][l] - phPad[k+1][l+1]);
          D2 = funcwrap1(phPad[k+1][l-1] - phPad[k][l]) -   funcwrap1(phPad[k][l] - phPad[k-1][l+1]);
          phaseQualityMap[k][l] = sqrt(H*H + V*V + D1*D1 + D2*D2);
          if(max<phaseQualityMap[k][l])max=phaseQualityMap[k][l];
	}
 }
	//resuing the given array 
	 for(int i=1;i<=m;i++){
		for(int j=1;j<=n;j++) phaseQualityMap[i-1][j-1]=phPad[i][j]/max;
	}
	
	for(int i=0;i<n;i++){
		ph[i][0]=1;
		ph[i][1]=1;
		ph[i][m-2]=1;
		ph[i][m-1]=1;
	}
	for(int i=0;i<m;i++){
		ph[0][i]=1;
		ph[1][i]=1;
		ph[n-1][i]=1;
		ph[n-2][i]=1;
	}
	return ph;
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
	vector<phase> sort_values; //sorted values of quality map; 
	 for(int i=0;i<n;i++){
	 	for(int j=0;j<m;j++){
	 		phase t;
	 		t.val = phaseQualityMap[i][j];
	 		t.row=i;
	 		t.col = j;
	 		sort_values.push_back(t);
		 }
	 }
	 sort(sort_values.begin(),sort_values.end(),comp);
	 vector<double> histcounts = linespace(0,1,nbins);
	 vector<vector<int> > treated (n,vector<int>(m,0)); //which pixels were treated; 
	 vector<vector<double> > phase2DQG (n,vector<double>(m,0));//final ans;
	 for(int i=0;i<histcounts.size();i++){
	 	double thr = histcounts[i];
	 	funcQGPU_Hist(ph,phaseQualityMap,treated,sort_values,thr,phase2DQG);
	 }
}
int main(){
	int m;
	int n;
	vector<vector<int> > a(n,vector<int>(m));
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			cin>>a[i][j];
		}
	}
	
	
	std::vector<std::vector<double>> x, y, z;
    for (double i = -5; i <= 5;  i += 0.25) {
        std::vector<double> x_row, y_row, z_row;
        for (double j = -5; j <= 5; j += 0.25) {
            x_row.push_back(i);
            y_row.push_back(j);
            z_row.push_back(::std::sin(::std::hypot(i, j)));
        }
        x.push_back(x_row);
        y.push_back(y_row);
        z.push_back(z_row);
    }

    plt::plot_surface(x, y, z);
    plt::show();
	
	
	
	
	
	
	return 0;
	
}
