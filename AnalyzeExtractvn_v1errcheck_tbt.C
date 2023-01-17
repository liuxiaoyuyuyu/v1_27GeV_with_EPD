//R__LOAD_LIBRARY(../../EPDsimulator/StRoot/StEpdUtil/StEpdGeom_cxx.so)

//#define export_txt
#define export_root //save the combined v1 TGraphErrors

using namespace std;

//void AnalyzeExtractvn_vzbinbybin(TString mOutputFileName="v1_16vz_Dec22_pol3odd",TString infilename=Form("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/fitoutput/root/Extracted_Nov25Fit_6f7ba7ae3e85c61517caa7c0fb2e82d336f24b68.root"),TString resolutionFile="/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/fitoutput/root/FlowHists_v1mike_run12.root"){
void AnalyzeExtractvn_v1errcheck_tbt(
//TString mOutputFileName="rbr_v1errcheck_pol5odd_run26_pretty",
TString mOutputFileName="tbt_v1errcheck_pol5odd_run26_extra_pretty_ver2_rainbow",
TString infilename="/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/fitoutput/root/Analysis_run26_v1_combinedTT_rmBadFits_demo2.root",
TString resolutionFile="/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/FlowHists_run26.root"){
//TString mOutputFileName="tbt_v1errcheck_pol5odd_run24",
//TString infilename=Form("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/fitoutput/root/Analysis_run24_v1_combinedTT_rmBadFits.root"),
//TString resolutionFile="/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/fitoutput/root/FlowHists_v1mike_run24.root"){
    double labSize(0.08),titSize(0.08);
    gStyle->SetPalette(kBlueGreenYellow);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleSize(labSize,"xy");
    gStyle->SetTitleSize(labSize,"t");
    gStyle->SetTitleAlign(23);  // first number is horizontal 1=right/2=center/3=left  second number is vertical
    gStyle->SetLabelSize(labSize,"xy");
    //gStyle->SetPadTopMargin(0.01);
    //gStyle->SetPadBottomMargin(0.15);
    //gStyle->SetPadLeftMargin(0.13);
    //gStyle->SetPadRightMargin(0.01);
    //gStyle->SetStatX(0.65);   // horizontal location of upper right corner  (in fraction of pad)
    //gStyle->SetStatY(0.4);   // vertical location of upper right corner  (in fraction of pad)
    //gStyle->SetStatW(0.2);   // width (but in what units???  not fraction of pad.  okay, it seems like fraction of pad / 2 or something...
    //gStyle->SetStatH(0.25);   // height (but in what units???  not fraction of pad
    gROOT->ForceStyle();
    TString OutputPDFName = "/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/pdf/"+mOutputFileName+".pdf";
    //StEpdGeom* geo = new StEpdGeom();

    TFile* tf = new TFile(infilename.Data(),"READ");
    TFile* ref = new TFile(resolutionFile.Data(),"READ");

    TProfile* Resolution[9][16];//cent, vz
    TProfile2D* AveEta[9][2];//cent, ew
    TGraphErrors* vnRaw[9][2][16];//cent, ew, vz 
    TGraphErrors* mgCom[9][16];//cent, vz

    double pi= TMath::Pi();
    double etabinsw[9]={2.1,2.4,2.8,3.1,3.4,3.7,4.0,4.3,5.1};
    double etabinse[9]={-5.1,-4.3,-4.0,-3.7,-3.4,-3.1,-2.8,-2.4,-2.1};

    for(int cent=0;cent<9;cent++){
        for(int ivz=0;ivz<16;ivz++){
            //Resolution[cent][ivz]  = (TProfile*)ref->Get(Form("Resolution%dPsi1Cent%dVz%d",0,cent,ivz));
            Resolution[cent][ivz]  = (TProfile*)ref->Get(Form("ResolutionPsi0Cent%dVz%d",cent,ivz));
            for(int ew=0;ew<2;ew++){
                //vnRaw[cent][ew][ivz] = (TGraphErrors*)tf->Get(Form("v1Cent%dEW%dVz%dEta%d",cent,ew,ivz,0));
                //vnRaw[cent][ew][ivz] = (TGraphErrors*)tf->Get(Form("v1Cent%dEW%dVz%dTT24",cent,ew,ivz));
                vnRaw[cent][ew][ivz] = (TGraphErrors*)tf->Get(Form("v1Cent%dEW%dVz%dTTAll",cent,ew,ivz));
            }
        }
        for(int ew=0;ew<2;ew++){
            AveEta[cent][ew] = (TProfile2D*)ref->Get(Form("AveEtaCent%dEW%d",cent,ew));
        }
    }

    //the following arrays are for making vn(EPD) graphs 
    double x[9][2][16][16];
    double y[9][2][16][16];
    double xerr[9][2][16][16];
    double yerr[9][2][16][16];//cent, order, ew, vz, ring
    
    //32 groups of points 
    double xCom[9][32][16];
    double yCom[9][32][16];
    double xerrCom[9][32][16]={0.0};
    double yerrCom[9][32][16];//cent, order, ew, vz, ring
    
    //16 vz bins 
    double xV[9][16][32];
    double yV[9][16][32];
    double xVerr[9][16][32];
    double yVerr[9][16][32];//cent, order, ew, vz, ring
    
    //16 vz bins together
    double xAll[9][512];
    double yAll[9][512];
    double xerrAll[9][512]={0.0};
    double yerrAll[9][512];//cent, order, ew, vz, ring
    
    //double color[16]={51,61,66,6,80,87,91,94,2,28,25,30,39,40,46,71}; 
    //int color[16]={51,861,66,6,829,819,91,94,2,28,25,30,39,863,46,433}; 
    int color[16]={634,632,910,896,807,801,800,829,819,419,839,432,433,861,860,881};
 
    TLine* l1=new TLine(-5.5,0.0,5.5,0.0);
    l1->SetLineStyle(2);
    TLine* l2=new TLine(-3.4,-0.13,-3.4,0.13);
    l2->SetLineStyle(2);
    TLine* l3=new TLine(3.4,-0.13,3.4,0.13);
    l3->SetLineStyle(2);
    TLine* l4=new TLine(0.0,-0.13,0.0,0.13);
    l4->SetLineStyle(2);
    TLine* l5=new TLine(2.0,0.0,5.5,0.0);
    l5->SetLineStyle(2);
    TLine* l6=new TLine(3.4,-0.019,3.4,0.019);
    l6->SetLineStyle(2);

    double edge[17] = {4.6, 9.0, 13.4, 17.8, 23.33, 28.86, 34.39, 39.92, 45.45, 50.98, 56.51, 62.05, 67.58, 73.11, 78.64, 84.17, 89.70};  // radii of ring boundaries
    double zEPD=375.0;  // epd distance to center of STAR
    double vzrange[17]={-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40};
    int vzrangeint[17]={-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40};
    int centrange[10]={80,70,60,50,40,30,20,10,5,0};
    
    double xCloud[9][32];//cent,32 groups
    double yCloud[9][32];  
    double xerrCloud[9][32]={0.0};  
    double yerrCloud[9][32];  

    for(int cent=0;cent<9;cent++){
        vector<tuple<double, double, double,double>> v;
        for(int ivz=0;ivz<16;ivz++){
            double res_TPC = TMath::Sqrt(Resolution[cent][ivz]->GetBinContent(1)*Resolution[cent][ivz]->GetBinContent(2)/Resolution[cent][ivz]->GetBinContent(3));
            for(int ew=0;ew<2;ew++){
                for(int iring=1;iring<=16;iring++){
                    //AveEta
                    x[cent][ew][ivz][iring-1]=AveEta[cent][ew]->GetBinContent(iring,ivz+1);
                    xerr[cent][ew][ivz][iring-1]=AveEta[cent][ew]->GetBinError(iring,ivz+1);
                    //xerr[cent][ew][ivz][iring-1]=0.0; 
                    
                    y[cent][ew][ivz][iring-1]=(vnRaw[cent][ew][ivz]->GetY())[iring-1];
                    y[cent][ew][ivz][iring-1]/=res_TPC;
                    yerr[cent][ew][ivz][iring-1]=vnRaw[cent][ew][ivz]->GetErrorY(iring-1)/res_TPC;
                    //IMPORTANT:GetErrorY(i), i starts at 0!!!!!!!!!ZERO ZERO ZERO!!!!!!!!
                    //yerr[cent][iorder][ew][ivz][iring-1]=0.0;
                    //cout<<"ring "<<iring<<": x="<<x[cent][iorder][ew][ivz][iring-1]<<", y="<<y[cent][iorder][ew][ivz][iring-1]<<", xerr="<<xerr[cent][iorder][ew][ivz][iring-1]<<", yerr="<<yerr[cent][iorder][ew][ivz][iring-1]<<endl;
                    v.push_back(make_tuple(x[cent][ew][ivz][iring-1],y[cent][ew][ivz][iring-1],xerr[cent][ew][ivz][iring-1],yerr[cent][ew][ivz][iring-1]));
                    //cout<<"x: "<<x[cent][iorder][ew][ivz][iring-1]<<endl;
                    xV[cent][ivz][ew*16+iring-1]=x[cent][ew][ivz][iring-1];
                    xVerr[cent][ivz][ew*16+iring-1]=xerr[cent][ew][ivz][iring-1];
                    yV[cent][ivz][ew*16+iring-1]=y[cent][ew][ivz][iring-1];
                    yVerr[cent][ivz][ew*16+iring-1]=yerr[cent][ew][ivz][iring-1];

                }//ring  

            }//ew
            
            mgCom[cent][ivz] = new TGraphErrors(32,xV[cent][ivz],yV[cent][ivz],xVerr[cent][ivz],yVerr[cent][ivz]);
            mgCom[cent][ivz]->SetMaximum(0.17);
            mgCom[cent][ivz]->SetMinimum(-0.17);
            mgCom[cent][ivz]->SetMarkerColor(color[ivz]);
            mgCom[cent][ivz]->SetLineColor(color[ivz]);
            mgCom[cent][ivz]->SetMarkerStyle(20);
            mgCom[cent][ivz]->SetMarkerSize(0.8);
            mgCom[cent][ivz]->SetTitle(Form("%d~%d%s",centrange[cent+1],centrange[cent],"%"));

        }//ivz
        sort(v.begin(),v.end());
        for(int ig=0;ig<32;ig++){
            for(int id=0;id<16;id++){
                xCom[cent][ig][id]=get<0>(v[ig*16+id]);
                yCom[cent][ig][id]=get<1>(v[ig*16+id]);
                xerrCom[cent][ig][id]=get<2>(v[ig*16+id]);
                yerrCom[cent][ig][id]=get<3>(v[ig*16+id]);
                //weird: after some centrality loops, yerrCom[0][0][0-1][0-15] becomes xCom
                //i.e eta values instead of the error bars. Therefore I am calculating yerrCloud right away
                
                xAll[cent][ig*16+id]=get<0>(v[ig*16+id]);
                yAll[cent][ig*16+id]=get<1>(v[ig*16+id]);
                xerrAll[cent][ig*16+id]=get<2>(v[ig*16+id]);
                yerrAll[cent][ig*16+id]=get<3>(v[ig*16+id]);
            }
            double xsum=0.0;
            double ysum=0.0;
            double err2sum=0.0;
            for(int id=0;id<16;id++){
                xsum+=xCom[cent][ig][id];
                ysum+=yCom[cent][ig][id]; 
                err2sum+=pow(yerrCom[cent][ig][id],2);
            }
            xsum/=16.0;
            ysum/=16.0;
            err2sum=sqrt(err2sum)/16.0;
            xCloud[cent][ig]=xsum;
            yCloud[cent][ig]=ysum;
            yerrCloud[cent][ig]=err2sum;
            //cout<<"cent="<<cent<<" ig="<<ig<<" yerrCloud="<<yerrCloud[cent][ig]<<endl;
            }
        v.clear();
    }//cent
  
    unique_ptr<TCanvas> TitlePage(new TCanvas("Title","Title",1400,1000));
    TPaveText tptitle(0.1,0.1,0.8,0.8);
    tptitle.AddText("Fitting all 16 vz bins and show the residuals, rebin the graph to 32 points and show the res.");
    tptitle.Draw();
    TitlePage->SaveAs(Form("%s%s",OutputPDFName.Data(),"("));//cannot .reset() for now because it will be saved again at the end
    
    TCanvas* cTot[2];//for plotting all the centralities 
    cTot[0]= new TCanvas("cTot0","cTot0",2000,1000);
    cTot[0]->Divide(4,2,0,0);
    cTot[1]= new TCanvas("cTot1","cTot1",2000,1000);
    cTot[1]->Divide(4,2,0,0);

    for(int cent=1;cent<9;cent++){
        cTot[0]->cd(cent);
        mgCom[cent][0]->Draw("AP");
        for(int ivz=1;ivz<16;ivz++){
            mgCom[cent][ivz]->Draw("P");
        }
    }

    
    //TH1D* frame[9];
    //TH1D* frame2[9];
    TCanvas* c1[9];//cent
    TCanvas* c2[9];//cent
    TCanvas* c3[9];//cent,different ways of combing points
    TCanvas* c4[9];
    TCanvas* c5[9];//All color
    /*
    double xCloud[9][32];//cent,32 groups
    double yCloud[9][32];  
    double xerrCloud[9][32]={0.0};  
    double yerrCloud[9][32];  
    */
    double xResCloud[9][32];  
    double yResCloud[9][32];  
    TGraphErrors* gvnCloud[9];
    TGraph* gvnRes[9];
    TGraphErrors* gvnAll[9];
    TGraph* gvnAllRes[9];
    TGraph* gvnAllvzRes[9][16];
    double xResAll[9][512];
    double yResAll[9][512];
    double xResAllvz[9][16][32];
    double yResAllvz[9][16][32];

    TH1D* ResAllProj[9];
    TH1D* ResCloudProj[9];
    
    double lex[4]={0.6,0.35,0.15,0.35}; 
    double ley[4]={0.05,0.05,0.6,0.6}; 
    double leh=0.25;
    double lew=0.1;
    TLegend* le[4];
    for(int i=0;i<4;i++){
        le[i]=new TLegend(lex[i],ley[i],lex[i]+lew,ley[i]+leh);
        for(int j=i*4;j<i*4+4;j++){
            le[i]->AddEntry(mgCom[0][j],Form("%d<vz<%d cm",vzrangeint[j],vzrangeint[j+1]),"PL");
            le[i]->SetBorderSize(0);
            le[i]->SetTextSize(0.045);
        }
    }

    for(int cent=0;cent<9;cent++){
        ResAllProj[cent] = new TH1D(Form("ResAllProjCent%d",cent),Form("ResAllProjCent%d",cent),20,-10,10);
        ResCloudProj[cent] = new TH1D(Form("ResCloudProjCent%d",cent),Form("ResCloudProjCent%d",cent),20,-10,10);

        TString title[3]={"small #eta range linear fit","average #eta and v_{1}, error propagation","average #eta and v_{1},RMS/#sqrt{N}"};
        c3[cent]=new TCanvas(Form("v1AllProjCent%d",cent),Form("v1ProjCent%d",cent),800,600);
        //gStyle->SetOptFit(0);
        c1[cent]=new TCanvas(Form("v1AllCent%d",cent),Form("v1Cent%d",cent),600,600);
        //c1[cent]->Update();
        TPad *pad3 = new TPad("pad3","pad3",0,0.33,1,1);
        TPad *pad4 = new TPad("pad4","pad4",0,0,1,0.33);
        pad3->SetBottomMargin(0.00001);
        pad3->SetBorderMode(0);
        //pad3->SetLogy();
        pad4->SetTopMargin(0.00001);
        pad4->SetBottomMargin(0.1);
        pad4->SetBorderMode(0);
        pad3->Draw();
        pad4->Draw();

        gvnAll[cent] = new TGraphErrors(512,xAll[cent],yAll[cent],xerrAll[cent],yerrAll[cent]);
        gvnAll[cent]->SetMarkerSize(0.8);
        gvnAll[cent]->SetMarkerStyle(20);
        gvnAll[cent]->SetMaximum(0.13);
        gvnAll[cent]->SetMinimum(-0.13);
        gvnAll[cent]->SetTitle(Form("%d~%d%s pol5_odd+p0 fit",centrange[cent+1],centrange[cent],"%"));
        //gvnAll[cent]->SetTitle(Form("%d~%d%s pol11 fit",centrange[cent+1],centrange[cent],"%"));
        gvnAll[cent]->GetXaxis()->SetTitle("#eta");
        gvnAll[cent]->GetYaxis()->SetTitle("v_{1}(#eta)");
        gvnAll[cent]->GetXaxis()->SetLimits(-5.5,5.5);
        //TF1* f1 = new TF1("f1","[6]+[0]*x+[1]*x^3+[2]*x^5+[3]*x^7+[4]*x^9+[5]*x^11",-6,6);
        TF1* f1 = new TF1("f1","[3]+[0]*x+[1]*x^3+[2]*x^5",-6,6);
        //TF1* f1 = new TF1("f1","[3]+[0]*x+[1]*x^3+[2]*x^5+[3]*x^7",-6,6);
        //TF1* f1 = new TF1("f1","[0]+[1]*x+[2]*x^3+[3]*x^5+[4]*x^7+[5]*x^9+[6]*x^11+[7]*x^13+[8]*x^2+[9]*x^4",-6,6);
        //TF1* f1 = new TF1("f1","pol11",-6,6);
        gvnAll[cent]->Fit("f1","FQ");
        TF1 *myfunc = gvnAll[cent]->GetFunction("f1");
        myfunc->SetLineColor(1);
        for(int i=0;i<512;i++){
            xResAll[cent][i]=xAll[cent][i];
            yResAll[cent][i]=(yAll[cent][i]-(*f1)(xResAll[cent][i]))/yerrAll[cent][i];
            ResAllProj[cent]->Fill(yResAll[cent][i]);
        }
        gvnAllRes[cent]=new TGraph(512,xResAll[cent],yResAll[cent]);
        gvnAllRes[cent]->SetMarkerSize(0.8);
        gvnAllRes[cent]->SetMarkerStyle(20);
        gvnAllRes[cent]->SetMaximum(10);
        gvnAllRes[cent]->SetMinimum(-10);
        gvnAllRes[cent]->GetXaxis()->SetTitle("#eta");
        gvnAllRes[cent]->GetYaxis()->SetTitle("Residual");
        gvnAllRes[cent]->SetTitle(Form("%d~%d%s residual",centrange[cent+1],centrange[cent],"%"));
        gvnAllRes[cent]->GetXaxis()->SetLimits(-5.5,5.5);
        TLine* lres1=new TLine(-5.5,1.0,5.5,1.0);
        TLine* lres2=new TLine(-5.5,-1.0,5.5,-1.0);
        TLine* lres3=new TLine(-5.5,5.0,5.5,5.0);
        TLine* lres4=new TLine(-5.5,-5.0,5.5,-5.0);    
        TLine* lres5=new TLine(-5.5,0.0,5.5,0.0);    
        lres1->SetLineStyle(2);
        lres2->SetLineStyle(2);
        lres3->SetLineStyle(2);
        lres4->SetLineStyle(2);
        pad3->cd();
        gvnAll[cent]->Draw("AP");
        pad4->cd(); 
        gvnAllRes[cent]->Draw("AP");
        lres1->Draw("same");
        lres2->Draw("same");
        lres3->Draw("same");
        lres4->Draw("same");
        lres5->Draw("same");
        //c1[cent]->SaveAs(OutputPDFName);
        
        c5[cent]=new TCanvas(Form("v1ColorAllCent%d",cent),Form("v1ColorCent%d",cent),600,600);
        //gStyle->SetOptFit(0);
        //c5[cent]->Update();
        for(int ivz=0;ivz<16;ivz++){
            for(int ig=0;ig<32;ig++){
                xResAllvz[cent][ivz][ig]=xV[cent][ivz][ig];
                yResAllvz[cent][ivz][ig]=(yV[cent][ivz][ig]-(*f1)(xResAllvz[cent][ivz][ig]))/yVerr[cent][ivz][ig];
            }
            gvnAllvzRes[cent][ivz]=new TGraph(32,xResAllvz[cent][ivz],yResAllvz[cent][ivz]);
            gvnAllvzRes[cent][ivz]->SetMarkerSize(0.8);
            gvnAllvzRes[cent][ivz]->SetMarkerStyle(20);
            gvnAllvzRes[cent][ivz]->SetMaximum(10);
            gvnAllvzRes[cent][ivz]->SetMinimum(-10);
            gvnAllvzRes[cent][ivz]->GetXaxis()->SetTitle("#eta");
            gvnAllvzRes[cent][ivz]->GetYaxis()->SetTitle("Residual");
            gvnAllvzRes[cent][ivz]->SetTitle(Form("%d~%d%s residual",centrange[cent+1],centrange[cent],"%"));
            gvnAllvzRes[cent][ivz]->GetXaxis()->SetLimits(-5.5,5.5);
            gvnAllvzRes[cent][ivz]->SetMarkerColor(color[ivz]);
            gvnAllvzRes[cent][ivz]->SetLineColor(color[ivz]);
        }
        TPad *pad5 = new TPad("pad5","pad5",0,0.33,1,1);
        TPad *pad6 = new TPad("pad6","pad6",0,0,1,0.33);
        pad5->SetBottomMargin(0.00001);
        pad5->SetBorderMode(0);
        //pad3->SetLogy();
        pad6->SetTopMargin(0.00001);
        pad6->SetBottomMargin(0.1);
        pad6->SetBorderMode(0);
        pad5->Draw();
        pad6->Draw();
        pad5->cd();
        gvnAll[cent]->Draw("AP");
        for(int ivz=0;ivz<16;ivz++){
            mgCom[cent][ivz]->Draw("P");
            for(int i=0;i<4;i++){
                le[i]->Draw("same");
            }
        }

        pad6->cd(); 
        gvnAllRes[cent]->Draw("AP");
        for(int ivz=0;ivz<16;ivz++){
            gvnAllvzRes[cent][ivz]->Draw("P");
        }
        lres1->Draw("same");
        lres2->Draw("same");
        lres3->Draw("same");
        lres4->Draw("same");
        lres5->Draw("same");
        c5[cent]->SaveAs(OutputPDFName);
        
        //gStyle->SetOptFit(1);
        c3[cent]->Update();
        c3[cent]->cd();
        ResAllProj[cent]->SetTitle("");
        ResAllProj[cent]->GetXaxis()->SetTitle("Residual");
        ResAllProj[cent]->Draw();
        ResAllProj[cent]->Fit("gaus","FQ");
        c3[cent]->SaveAs(OutputPDFName);

        c4[cent]=new TCanvas(Form("v1CloudProjCent%d",cent),Form("v1ProjCent%d",cent),800,600);
        c2[cent]=new TCanvas(Form("v1CloudCent%d",cent),Form("v1Cent%d",cent),600,600);
        //gStyle->SetOptFit(0);
        //c2[cent]->Update();
        TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
        TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
        pad1->SetBottomMargin(0.00001);
        pad1->SetBorderMode(0);
        //pad1->SetLogy();
        pad2->SetTopMargin(0.00001);
        pad2->SetBottomMargin(0.1);
        pad2->SetBorderMode(0);
        pad1->Draw();
        pad2->Draw();
        /*
        for(int ig=0;ig<32;ig++){
            double xsum=0.0;
            double ysum=0.0;
            double err2sum=0.0;
            for(int id=0;id<16;id++){
               xsum+=xCom[cent][0][ig][id];
               ysum+=yCom[cent][0][ig][id]; 
               err2sum+=pow(yerrCom[cent][0][ig][id],2);
            }
            xsum/=16.0;
            ysum/=16.0;
            err2sum=sqrt(err2sum)/16.0;
            xCloud[cent][ig]=xsum;
            yCloud[cent][ig]=ysum;
            yerrCloud[cent][ig]=err2sum;
            //cout<<"cent="<<cent<<" ig="<<ig<<" yerrCloud="<<yerrCloud[cent][ig]<<endl;
        }
        */

        //add loop for different points combination 
        gvnCloud[cent] = new TGraphErrors(32,xCloud[cent],yCloud[cent],xerrCloud[cent],yerrCloud[cent]);
        gvnCloud[cent]->SetMarkerSize(0.8);
        gvnCloud[cent]->SetMarkerStyle(20);
        gvnCloud[cent]->SetMaximum(0.13);
        gvnCloud[cent]->SetMinimum(-0.13);
        gvnCloud[cent]->SetTitle(Form("%d~%d%s pol5_odd+p0 fit",centrange[cent+1],centrange[cent],"%"));
        //gvnCloud[cent]->SetTitle(Form("%d~%d%s pol11 fit",centrange[cent+1],centrange[cent],"%"));
        gvnCloud[cent]->GetXaxis()->SetTitle("#eta");
        gvnCloud[cent]->GetYaxis()->SetTitle("v_{1}(#eta)");
        gvnCloud[cent]->GetXaxis()->SetLimits(-5.5,5.5);
        
        //TF1* f2 = new TF1("f2","[3]+[0]*x+[1]*x^3+[2]*x^5+[3]*x^7",-6,6);
        //TF1* f2 = new TF1("f2","[6]+[0]*x+[1]*x^3+[2]*x^5+[3]*x^7+[4]*x^9+[5]*x^11",-6,6);
        //TF1* f2 = new TF1("f2","[0]*x+[1]*x^3+[2]*x^5+[3]*x^7+[4]*x^9+[5]*x^11",-6,6);
        TF1* f2 = new TF1("f2","[3]+[0]*x+[1]*x^3+[2]*x^5",-6,6);
        //TF1* f2 = new TF1("f2","pol11",-6,6);
        //gvnCloud[cent]->Fit("f2","FQ0");
        gvnCloud[cent]->Fit("f2","FQ");
        for(int i=0;i<32;i++){
            xResCloud[cent][i]=xCloud[cent][i];
            yResCloud[cent][i]=(yCloud[cent][i]-(*f2)(xResCloud[cent][i]))/yerrCloud[cent][i];
            ResCloudProj[cent]->Fill(yResCloud[cent][i]);
        }
        gvnRes[cent]=new TGraph(32,xResCloud[cent],yResCloud[cent]);
        gvnRes[cent]->SetMarkerSize(0.8);
        gvnRes[cent]->SetMarkerStyle(20);
        gvnRes[cent]->SetMaximum(10);
        gvnRes[cent]->SetMinimum(-10);
        gvnRes[cent]->SetTitle(Form("%d~%d%s residual",centrange[cent+1],centrange[cent],"%"));
        gvnRes[cent]->GetXaxis()->SetTitle("#eta");
        gvnRes[cent]->GetYaxis()->SetTitle("Residual");
        gvnRes[cent]->GetXaxis()->SetLimits(-5.5,5.5);
        
        pad1->cd();
        gvnCloud[cent]->Draw("AP");
        /*
        for(int ivz=0;ivz<16;ivz++){
            mgCom[cent][ivz]->Draw("P");
        }
        myfunc->Draw("same");
        gvnCloud[cent]->Draw("P");
        */
        pad2->cd(); 
        gvnRes[cent]->Draw("AP");
        /*
        for(int ivz=0;ivz<16;ivz++){
            gvnAllvzRes[cent][ivz]->Draw("P");
        }
        gvnRes[cent]->Draw("P");
        */
        lres1->Draw("same");
        lres2->Draw("same");
        lres3->Draw("same");
        lres4->Draw("same");
        lres5->Draw("same");
        c2[cent]->SaveAs(OutputPDFName);
        //gStyle->SetOptFit(1);
        c4[cent]->Update();
        c4[cent]->cd();
        ResCloudProj[cent]->SetTitle("");
        ResCloudProj[cent]->GetXaxis()->SetTitle("Residual");
        ResCloudProj[cent]->Draw();
        ResCloudProj[cent]->Fit("gaus","FQ");
        c4[cent]->SaveAs(OutputPDFName);

    }
    
    //gStyle->SetOptFit(0);
    //cTot[1]->Update();
    
    for(int cent=1;cent<9;cent++){
        cTot[1]->cd(cent);
        TF1 *myfunc = gvnCloud[cent]->GetFunction("f2");
        myfunc->Draw("0");
        gvnCloud[cent]->Draw("AP");
        
    }
    for(int i=0;i<2;i++){
        cTot[i]->SaveAs(OutputPDFName);
    }

    TitlePage->SaveAs(Form("%s%s",OutputPDFName.Data(),"]"));//cannot .reset() for now because it will be saved again at the end

    #ifdef export_root
    TFile* froot= new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/root/tbt_v1_run26_combineVz_ver2.root","RECREATE");
    froot->cd();
    for(int icent=0;icent<9;icent++){
        gvnCloud[icent]->SetName(Form("v1Cent%d",icent));
        gvnCloud[icent]->Write();
    }
    froot->Write();
    froot->Close();
    #endif

    #ifdef export_txt
    std::ofstream flowlist;
    flowlist.open(Form("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/datapoints/FlowDataPoints_%s.txt",mOutputFileName.Data()));   
    std::string contents;
    flowlist<<"v1(EPD) at AuAu 27GeV, nMip fitting (Extracted_Dec23_v1_onlycov.root), errs from the covariance matrix."<<endl;
    flowlist<<"double xv1EPDEW[9][32]={";
    for(int cent=0;cent<9;cent++){
        for(int i=0;i<32;i++){
            contents.append(to_string(xCloud[cent][i]));
            contents.append(",");
        }
    }
    contents.pop_back();
    flowlist<<contents<<"};"<<endl;
    contents.clear();
    flowlist<<"double yv1EPDEW[9][32]={";
    for(int cent=0;cent<9;cent++){
        for(int i=0;i<32;i++){
            contents.append(to_string(yCloud[cent][i]));
            contents.append(",");
        }
    }
    contents.pop_back();
    flowlist<<contents<<"};"<<endl;
    contents.clear();
    flowlist<<"double xerrv1EPDEW[9][32]={";
    for(int cent=0;cent<9;cent++){
        for(int i=0;i<32;i++){
            contents.append(to_string(xerrCloud[cent][i]));
            contents.append(",");
        }
    }
    contents.pop_back();
    flowlist<<contents<<"};"<<endl;
    contents.clear();
    flowlist<<"double yerrv1EPDEW[9][32]={";
    for(int cent=0;cent<9;cent++){
        for(int i=0;i<32;i++){
            contents.append(to_string(yerrCloud[cent][i]));
            contents.append(",");
        }
    }
    contents.pop_back();
    flowlist<<contents<<"};"<<endl;
    contents.clear();
    flowlist.close();
    #endif

}
