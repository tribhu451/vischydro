// root -b -q "centrality_determin.C(\"au1.root\")"

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"

using namespace std;

void centrality_determin(const char* aname)
 {
  
   
   double nch; double imp_b; int npart;
   
   
   // inputs [
   int bins=8000;
   double min_bin_val = 0.;
   double max_bin_val = 810.;
   
   double xcut[10]={0.0};
   
   const int pr_f = 7;
   double percent[pr_f] = {0,  5,  15, 25 , 35, 45, 55};
   
   TFile* f1=new TFile(aname,"read");
   TH1F* h1 = new TH1F("h1","",bins,min_bin_val,max_bin_val);
   TTree* TreeAu =new TTree();
   f1->GetObject("Tree_opt_au",TreeAu);
   TreeAu->SetBranchAddress("NCh",&nch);
   TreeAu->SetBranchAddress("Impact_parameter",&imp_b);
   TreeAu->SetBranchAddress("N_Part",&npart);
   // inputs ]
   
   
   
   
   
   double maxch = 0 ; 
   int Entries=TreeAu->GetEntries();        
   for(int j=0;j<Entries;j++)
     {
       TreeAu->GetEntry(j);
       if(nch > max_bin_val)
	 {cout<<"nch > maximum bin range : "<<nch<<endl;}
       h1->Fill(nch);
       if (nch > maxch) maxch = nch ; 
     }
   
   int n1=0;
   double Sum=0.0;
   int per_count = 0 ;
   xcut[n1] = maxch ;
   
   for(int j=Entries;j>=1; j--)
     {
       double aSum = h1->GetBinContent(j);
       double perc=percent[per_count+1] - percent[per_count];   // perc=10 means 0-10% , 10-20%, 20-30%....                                     
       if(Sum>(Entries*perc*0.01))
	 {
           n1=n1+1; 
	   double p0=h1->GetBinCenter(j);                                       
	   xcut[n1]=p0;
           Sum=0.0;
	   
           per_count += 1;

           if(n1==pr_f){break;}  
         }
       else
	 {
	   Sum=Sum+aSum;
	 }
     }
   
   //p0=h1->GetBinCenter(1);                 
   //n1=n1+1;
   //xcut[n1]=p0;
   
   for(int i=0; i<pr_f-1; i++)
     {
       cout<<percent[i]<<"-"<<percent[i+1]<<" \%"<<endl;
     }
   
   for(int i=0 ; i<n1; i++)
     {
       if (i!=0)
         cout<<xcut[i-1]<<" - "<<xcut[i]<<endl;
     }
   
   
   double bb[n1];
   int npartt[n1];
   for(int j=0;j<Entries;j++)
     {
       TreeAu->GetEntry(j);
       for (int i = 0; i < n1 ; i ++){
	 if(nch < xcut[i]+5 && nch > xcut[i]-5)
	   { bb[i] = imp_b; npartt[i] = npart; }
	 
       }
     }
   
   for (int i = 0; i < n1 ; i ++)
     {
       cout<<"multiplicity "<<xcut[i]<<"\t=>\tb : "<<bb[i]<<"\tnpart : "<<npartt[i]<<endl;
     }
   
   
   
 }
