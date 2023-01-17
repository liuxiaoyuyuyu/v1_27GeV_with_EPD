//https://root.cern/doc/v610/canvas2_8C.html
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;
   for (Int_t i=0;i<Nx;i++) {
      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
      for (Int_t j=0;j<Ny;j++) {
         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            //vmard = 0.0;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            //vmaru = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
         C->cd(0);
         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
         pad->Draw();
      }
   }
}

//void DrawPrettyPlots_new(TString mOutputFileName="DrawPrettyPlots_new_run26"){
void DrawPrettyPlots_new_ver2(TString mOutputFileName="DrawPrettyPlots_new_run26_tbt_ver2"){
//void DrawPrettyPlots_new(TString mOutputFileName="DrawPrettyPlots_new_run26_rbr_GEANTcorrected"){
//void DrawPrettyPlots_new(TString mOutputFileName="DrawPrettyPlots_new_nonflow"){
    gStyle->SetTitleSize(0.08,"t");
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();
    //TFile* fStat = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/GEANT2021/root/DrawCorrectedv1_from_run9_03102022.root","READ");
    TFile* fStat = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/root/tbt_v1_run26_combineVz_ver2.root","READ");
    //TFile* fStat = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/root/rbr_v1_run26_combineVz.root","READ");
    //TFile* fStat = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/root/DrawCorrectedv1_RCFrun26_default_Sep102022.root","READ");
    //TFile* fStat = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/fitoutput/root/rbr_v1_nonflow_combineVz.root","READ");
    //TFile* fStat = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/fitoutput/root/rbr_v1_nonflow_combineVz_TPC0p8.root","READ");
    //TFile* fSys = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/GEANT2021/root/DrawSysCheck_add180rot.root","READ");
    //TFile* fSys = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/root/tbt_v1_run26_combineVz.root","READ");
    //TFile* fUrQMD = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_code/UrQMD/UrQMDAna/UrQMD_v1_allEP.root","READ");

    //TString OutputPDFName = "/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/GEANT2021/pdf/"+mOutputFileName+".pdf";
    //TString OutputPDFName = "/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/pdf/"+mOutputFileName+".pdf";
    TString OutputPDFName = "/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/fitoutput/pdf/"+mOutputFileName+".pdf";
    unique_ptr<TCanvas> TitlePage(new TCanvas("Title","Title",450,1200));
    TPaveText tptitle(0.1,0.1,0.8,0.8);
    tptitle.AddText(Form("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/DrawPrettyPlots_new_ver2.C"));
    tptitle.Draw();
    TitlePage->SaveAs(Form("%s%s",OutputPDFName.Data(),"("));//cannot .reset() for now because it will be saved again at the end

    TGraphErrors* Gv1Stat[9];
    //TGraphErrors* Gv1Sys[9];
    //TGraphErrors* Gv1UrQMD[9][2];//cent,ew
    //TProfile* GUrQMDcosphi[9];

    double mMarkerSize=0.9;
    for(int icent=0;icent<9;icent++){
        /*
        if (icent==0||icent==1){
            Gv1Stat[icent]=(TGraphErrors*)fStat->Get(Form("v1CorrectedCh1Cent%d",2));
        }
        */
        //Gv1Stat[icent]=(TGraphErrors*)fStat->Get(Form("v1CorrectedCh1Cent%d",icent));
        //Gv1Sys[icent]=(TGraphErrors*)fSys->Get(Form("v1CorrectedSysErrorCh1Cent%d",icent));
        Gv1Stat[icent]=(TGraphErrors*)fStat->Get(Form("v1Cent%d",icent));
        //Gv1Sys[icent]=(TGraphErrors*)fSys->Get(Form("v1Cent%d",icent));

        /*
        GUrQMDcosphi[icent]=(TProfile*)fUrQMD->Get(Form("CosphiCent%d",icent));
        GUrQMDcosphi[icent]->SetLineColor(2);
        GUrQMDcosphi[icent]->SetMarkerColor(2);
        GUrQMDcosphi[icent]->SetFillColor(2);
        GUrQMDcosphi[icent]->SetFillStyle(1001);
        for(int ew=0;ew<2;ew++){
            Gv1UrQMD[icent][ew]=(TGraphErrors*)fUrQMD->Get(Form("v1EPCent%dEW%dITPC13",icent,ew));
            Gv1UrQMD[icent][ew]->SetFillColorAlpha(601,0.5);
            Gv1UrQMD[icent][ew]->SetLineColorAlpha(601,1);
            Gv1UrQMD[icent][ew]->SetFillStyle(1001);
        }
        */
    }

    //Combine centralities. 0:60~80% (0+1); 1: 40~60% (2+3); 2: 20~40% (4+5); 3: 0~20% (0.5*7+0.5*8+6)
    double ax[4][2][16];
    double ay[4][2][16];
    double aexstat[4][2][16];
    //double aexsys[4][2][16];
    double aeystat[4][2][16];
    //double aeysys[4][2][16];
    
    double axMirror[4][2][16];
    double ayMirror[4][2][16];

    /*
    double xurqmd[3][2][16];
    double yurqmd[3][2][16];
    double xerrurqmd[3][2][16];
    double yerrurqmd[3][2][16];

    double xurqmdMirror[3][2][16];
    double yurqmdMirror[3][2][16];

    double xcos[3][160];
    double ycos[3][160];
    double xerrcos[3][160];
    double yerrcos[3][160];
    */

    for(int ew=0;ew<2;ew++){
        for(int ieta=0;ieta<16;ieta++){
            for(int icent=0;icent<4;icent++){
                ax[icent][ew][ieta]=Gv1Stat[2]->GetX()[16*ew+ieta];
                aexstat[icent][ew][ieta]=0.0;
                //aexsys[icent][ew][ieta]=0.05;
                //xurqmd[icent][ew][ieta]=Gv1UrQMD[2][ew]->GetX()[ieta];
                //xerrurqmd[icent][ew][ieta]=0.0; 
            }
            ay[0][ew][ieta]=(Gv1Stat[0]->GetY()[16*ew+ieta]+Gv1Stat[1]->GetY()[16*ew+ieta])/2.0;
            ay[1][ew][ieta]=(Gv1Stat[2]->GetY()[16*ew+ieta]+Gv1Stat[3]->GetY()[16*ew+ieta])/2.0;
            ay[2][ew][ieta]=(Gv1Stat[4]->GetY()[16*ew+ieta]+Gv1Stat[5]->GetY()[16*ew+ieta])/2.0;
            ay[3][ew][ieta]=(Gv1Stat[7]->GetY()[16*ew+ieta]+Gv1Stat[8]->GetY()[16*ew+ieta]+2.0*Gv1Stat[6]->GetY()[16*ew+ieta])/4.0;
        
            aeystat[0][ew][ieta]=0.5*sqrt(pow(Gv1Stat[0]->GetErrorY(16*ew+ieta),2)+pow(Gv1Stat[1]->GetErrorY(16*ew+ieta),2));
            aeystat[1][ew][ieta]=0.5*sqrt(pow(Gv1Stat[2]->GetErrorY(16*ew+ieta),2)+pow(Gv1Stat[3]->GetErrorY(16*ew+ieta),2));
            aeystat[2][ew][ieta]=0.5*sqrt(pow(Gv1Stat[4]->GetErrorY(16*ew+ieta),2)+pow(Gv1Stat[5]->GetErrorY(16*ew+ieta),2));
            aeystat[3][ew][ieta]=0.25*sqrt(pow(Gv1Stat[7]->GetErrorY(16*ew+ieta),2)+pow(Gv1Stat[8]->GetErrorY(16*ew+ieta),2)
            +4.0*pow(Gv1Stat[6]->GetErrorY(16*ew+ieta),2));
            
            /*
            aeysys[0][ew][ieta]=0.5*sqrt(pow(Gv1Sys[2]->GetErrorY(16*ew+ieta),2)+pow(Gv1Stat[3]->GetErrorY(16*ew+ieta),2));
            aeysys[1][ew][ieta]=0.5*sqrt(pow(Gv1Sys[4]->GetErrorY(16*ew+ieta),2)+pow(Gv1Stat[5]->GetErrorY(16*ew+ieta),2));
            aeysys[2][ew][ieta]=0.25*sqrt(pow(Gv1Sys[7]->GetErrorY(16*ew+ieta),2)+pow(Gv1Stat[8]->GetErrorY(16*ew+ieta),2)
            +4.0*pow(Gv1Sys[6]->GetErrorY(16*ew+ieta),2));

            yurqmd[0][ew][ieta]=(Gv1UrQMD[2][ew]->GetY()[ieta]+Gv1UrQMD[3][ew]->GetY()[ieta])/2.0;
            yurqmd[1][ew][ieta]=(Gv1UrQMD[4][ew]->GetY()[ieta]+Gv1UrQMD[5][ew]->GetY()[ieta])/2.0;
            yurqmd[2][ew][ieta]=(Gv1UrQMD[7][ew]->GetY()[ieta]+Gv1UrQMD[8][ew]->GetY()[ieta]+2.0*Gv1UrQMD[6][ew]->GetY()[ieta])/4.0;
            
            yerrurqmd[0][ew][ieta]=0.5*sqrt(pow(Gv1UrQMD[2][ew]->GetErrorY(ieta),2)+pow(Gv1UrQMD[3][ew]->GetErrorY(ieta),2));
            yerrurqmd[1][ew][ieta]=0.5*sqrt(pow(Gv1UrQMD[4][ew]->GetErrorY(ieta),2)+pow(Gv1UrQMD[5][ew]->GetErrorY(ieta),2));
            yerrurqmd[2][ew][ieta]=0.25*sqrt(pow(Gv1UrQMD[7][ew]->GetErrorY(ieta),2)+pow(Gv1UrQMD[8][ew]->GetErrorY(ieta),2)
            +4.0*pow(Gv1UrQMD[6][ew]->GetErrorY(ieta),2));
            */
        }
    }

    //TGraphErrors* geSysMirror[3][2];
    TGraphErrors* geStatMirror[4][2];
    //TGraphErrors* geUrQMDMirror[3][2];
    //Make the mirrored points 
    for(int icent=0;icent<4;icent++){
    //for(int icent=1;icent<4;icent++){
        for(int ew=0;ew<2;ew++){
            for(int ieta=0;ieta<16;ieta++){
                axMirror[icent][ew][ieta]=-ax[icent][ew][ieta];
                ayMirror[icent][ew][ieta]=-ay[icent][ew][ieta];
                //xurqmdMirror[icent][ew][ieta]=-xurqmd[icent][ew][ieta];
                //yurqmdMirror[icent][ew][ieta]=-yurqmd[icent][ew][ieta];
            }
            //geSysMirror[icent][ew]=new TGraphErrors(16,axMirror[icent][ew],ayMirror[icent][ew],aexsys[icent][ew],aeysys[icent][ew]);
            geStatMirror[icent][ew]=new TGraphErrors(16,axMirror[icent][ew],ayMirror[icent][ew],aexstat[icent][ew],aeystat[icent][ew]);
            geStatMirror[icent][ew]->SetMarkerStyle(24);
            geStatMirror[icent][ew]->SetMarkerSize(mMarkerSize);
            geStatMirror[icent][ew]->SetLineColor(2);
            geStatMirror[icent][ew]->SetMarkerColor(2);
            /*
            geSysMirror[icent][ew]->SetMarkerStyle(24);
            geSysMirror[icent][ew]->SetMarkerSize(mMarkerSize);
            geSysMirror[icent][ew]->SetLineColor(2);
            geSysMirror[icent][ew]->SetMarkerColor(2);
            geUrQMDMirror[icent][ew]=new TGraphErrors(16,xurqmdMirror[icent][ew],yurqmdMirror[icent][ew],xerrurqmd[icent][ew],yerrurqmd[icent][ew]);
            geUrQMDMirror[icent][ew]->SetFillColorAlpha(910,0.5);
            geUrQMDMirror[icent][ew]->SetLineColorAlpha(910,1);
            geUrQMDMirror[icent][ew]->SetFillStyle(1001);
            */
        }
    }         
    /*
    //combine UrQMD cosphi 
    for(int ieta=0;ieta<160;ieta++){
        for(int icent=0;icent<3;icent++){
            xcos[icent][ieta]=GUrQMDcosphi[2]->GetBinCenter(ieta+1);
            xerrcos[icent][ieta]=0.0;
        }
        ycos[0][ieta]=0.5*(GUrQMDcosphi[2]->GetBinContent(ieta+1)+GUrQMDcosphi[3]->GetBinContent(ieta+1));
        ycos[1][ieta]=0.5*(GUrQMDcosphi[4]->GetBinContent(ieta+1)+GUrQMDcosphi[5]->GetBinContent(ieta+1));
        ycos[2][ieta]=0.25*(GUrQMDcosphi[7]->GetBinContent(ieta+1)+GUrQMDcosphi[8]->GetBinContent(ieta+1)+2.0*GUrQMDcosphi[6]->GetBinContent(ieta+1));
        yerrcos[0][ieta]=0.5*sqrt(pow(GUrQMDcosphi[2]->GetBinError(ieta+1),2)+pow(GUrQMDcosphi[3]->GetBinError(ieta+1),2)); 
        yerrcos[1][ieta]=0.5*sqrt(pow(GUrQMDcosphi[4]->GetBinError(ieta+1),2)+pow(GUrQMDcosphi[5]->GetBinError(ieta+1),2)); 
        yerrcos[2][ieta]=0.25*sqrt(pow(GUrQMDcosphi[7]->GetBinError(ieta+1),2)+pow(GUrQMDcosphi[8]->GetBinError(ieta+1),2)
        +4.0*pow(GUrQMDcosphi[6]->GetBinError(ieta+1),2)); 
    }
    */
    int centint[5]={80,60,40,20,0};
    double metal[2]={-4.9,1.9};
    double metah[2]={-1.9,4.9};
    TGraphErrors* geStat[4][2];
    /*
    TGraphErrors* geSys[3][2];
    TGraphErrors* gUrQMD[3][2];
    TGraphErrors* gcos[3];
    */
    for(int icent=0;icent<4;icent++){
    //for(int icent=1;icent<4;icent++){
        for(int ew=0;ew<2;ew++){
            /*
            gUrQMD[icent][ew]=new TGraphErrors(16,xurqmd[icent][ew],yurqmd[icent][ew],xerrurqmd[icent][ew],yerrurqmd[icent][ew]);
            gUrQMD[icent][ew]->SetFillColorAlpha(601,0.5);
            gUrQMD[icent][ew]->SetLineColorAlpha(601,1);
            gUrQMD[icent][ew]->SetFillStyle(1001);
            */
            geStat[icent][ew]= new TGraphErrors(16,ax[icent][ew],ay[icent][ew],aexstat[icent][ew],aeystat[icent][ew]);
            geStat[icent][ew]->SetMarkerStyle(20);
            geStat[icent][ew]->SetMarkerSize(mMarkerSize);
            geStat[icent][ew]->SetLineColor(1);
            geStat[icent][ew]->SetMarkerColor(1);
            /*
            geSys[icent][ew]= new TGraphErrors(16,ax[icent][ew],ay[icent][ew],aexsys[icent][ew],aeysys[icent][ew]);
            geSys[icent][ew]->SetMarkerStyle(20);
            geSys[icent][ew]->SetMarkerSize(mMarkerSize);
            geSys[icent][ew]->SetLineColor(1);
            geSys[icent][ew]->SetMarkerColor(1);
            //geSys[icent][ew]->SetFillColorAlpha(1,1);
            geSys[icent][ew]->SetFillStyle(1001);
            */
            geStat[icent][ew]->SetMaximum(0.08);
            geStat[icent][ew]->SetMinimum(-0.17);
            geStat[icent][ew]->GetXaxis()->SetLimits(metal[ew],metah[ew]);
            geStat[icent][ew]->GetXaxis()->SetTitle("#eta");
            geStat[icent][ew]->GetXaxis()->SetTitleSize(0.08);
            geStat[icent][ew]->GetXaxis()->SetLabelSize(0.08);
            geStat[icent][ew]->GetXaxis()->SetTitleOffset(1);
            geStat[icent][ew]->GetXaxis()->CenterTitle();
            geStat[icent][ew]->GetYaxis()->SetTitle("v_{1}");
            geStat[icent][ew]->GetYaxis()->SetTitleSize(0.08);
            geStat[icent][ew]->GetYaxis()->SetLabelSize(0.08);
            geStat[icent][ew]->GetYaxis()->SetTitleOffset(0.8);
            //geStat[icent][ew]->GetYaxis()->CenterTitle();
            geStat[icent][ew]->SetTitle(Form("%d~%d%s",centint[icent+1],centint[icent],"%"));
            /*
            geSys[icent][ew]->SetMaximum(0.17);
            geSys[icent][ew]->SetMinimum(-0.17);
            geSys[icent][ew]->GetXaxis()->SetTitle("#eta");
            geSys[icent][ew]->GetXaxis()->CenterTitle();
            geSys[icent][ew]->GetYaxis()->SetTitle("v_{1}");
            geSys[icent][ew]->GetYaxis()->CenterTitle();
            geSys[icent][ew]->SetTitle(Form("%d~%d%s",centint[icent+1],centint[icent],"%"));
            */
        }
        /*
        gcos[icent]= new TGraphErrors(160,xcos[icent],ycos[icent],xerrcos[icent],yerrcos[icent]);
        gcos[icent]->SetLineColorAlpha(910,1);
        gcos[icent]->SetMarkerColorAlpha(910,1);
        gcos[icent]->SetFillColorAlpha(910,0.5);
        gcos[icent]->SetFillStyle(1001);
        */
    }

    double mybeam[2]={-3.4,3.4};
    TLine* l1[2];
    TLine* l2[2];
    for(int ew=0;ew<2;ew++){
        l1[ew]=new TLine(metal[ew],0.0,metah[ew],0.0);
        l2[ew]=new TLine(mybeam[ew],-0.17,mybeam[ew],0.17);
        l1[ew]->SetLineStyle(2);
        l2[ew]->SetLineStyle(10);
        l2[ew]->SetLineColorAlpha(801,1);
        l2[ew]->SetLineWidth(3);
    }
    TLine* l3 =new TLine(-4.9,0.0,4.9,0.0); 
    TLine* l4 =new TLine(0.0,-0.17,0.0,0.17);
    l3->SetLineStyle(2); 
    l4->SetLineStyle(2); 

    TLegend* le=new TLegend(0.3,0.05,0.9,0.25); 
    le->AddEntry(geStat[0][1],Form("#eta < 0"));
    le->AddEntry(geStatMirror[0][0],Form("#eta > 0 rotated 180 degrees"));
    le->SetBorderSize(0);
    le->SetTextSize(0.08);
    
    /*
    TLegend* le2=new TLegend(0.35,0.05,0.55,0.4); 
    le2->AddEntry(geStat[0][0],Form("STAR"));
    le2->AddEntry(gUrQMD[0][0],Form("UrQMD <cos(#phi-#Psi_{1}^{EP})>"));
    le2->AddEntry(gcos[0],"UrQMD <cos(#phi-#Psi_{1}^{RP})>","PL");
    le2->SetBorderSize(0);
    le2->SetTextSize(0.075);
    
    TCanvas *c2 = new TCanvas("c2", "c2",800,1200);
    c2->Range(-1.156662,-0.1237805,10.52709,1.103354);
    c2->SetFillColor(0);
    c2->SetBorderMode(0);
    c2->SetBorderSize(2);
    c2->SetRightMargin(0.01);
    c2->SetTopMargin(0.01);
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.15);
    c2->SetFrameBorderMode(0);
    c2->SetFrameBorderMode(0);
    c2->Divide(2,3,0,0);
    for(int icent=0;icent<3;icent++){
        for(int ew=0;ew<2;ew++){
            c2->cd(icent*2+ew+1);
            geStat[icent][ew]->Draw("AP");
            gUrQMD[icent][ew]->Draw("l3");
            //geSys[icent][ew]->SetFillColorAlpha(1,0.5);
            //geSys[icent][ew]->Draw("3");
            geSys[icent][ew]->Draw("PS;;5 s=0.5");
            geStat[icent][ew]->Draw("P");
            l1[ew]->Draw("same");
            l2[ew]->Draw("same");
        }
    }
    c2->SaveAs(OutputPDFName); 
    */
    /*
    // Margins
    Float_t lMargin = 0.12;
    Float_t rMargin = 0.05;
    Float_t bMargin = 0.15;
    Float_t tMargin = 0.05;

    TCanvas *C = (TCanvas*) gROOT->FindObject("C");
    if (C) delete C;
    C = new TCanvas("C","canvas",800,1200);
    C->SetFillStyle(4000);
    // Number of PADS
    const Int_t Nx = 2;
    const Int_t Ny = 3;
    // Canvas setup
    CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
    TPad *pad[Nx][Ny];
    for (Int_t i=0;i<Nx;i++) {
        for (Int_t j=0;j<Ny;j++) {
            C->cd(0);
            // Get the pads previously created.
            char pname[16];
            sprintf(pname,"pad_%i_%i",i,j);
            pad[i][j] = (TPad*) gROOT->FindObject(pname);
            pad[i][j]->Draw();
            pad[i][j]->SetFillStyle(4000);
            pad[i][j]->SetFrameFillStyle(4000);
            pad[i][j]->cd();
            // Size factors
            Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
            Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();
            // Format for y axis
            geStat[2-j][i]->GetYaxis()->SetLabelFont(43);
            geStat[2-j][i]->GetYaxis()->SetLabelSize(24);
            geStat[2-j][i]->GetYaxis()->SetLabelOffset(0.02);
            geStat[2-j][i]->GetYaxis()->SetTitleFont(43);
            geStat[2-j][i]->GetYaxis()->SetTitleSize(24);
            geStat[2-j][i]->GetYaxis()->SetTitleOffset(2);
            //geStat[2-j][i]->GetYaxis()->CenterTitle();
            geStat[2-j][i]->GetYaxis()->SetNdivisions(505);
            // TICKS Y Axis
            geStat[2-j][i]->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
            // Format for x axis
            geStat[2-j][i]->GetXaxis()->SetLabelFont(43);
            geStat[2-j][i]->GetXaxis()->SetLabelSize(24);
            geStat[2-j][i]->GetXaxis()->SetLabelOffset(0.02);
            geStat[2-j][i]->GetXaxis()->SetTitleFont(43);
            geStat[2-j][i]->GetXaxis()->SetTitleSize(24);
            geStat[2-j][i]->GetXaxis()->SetTitleOffset(2);
            geStat[2-j][i]->GetXaxis()->CenterTitle();
            geStat[2-j][i]->GetXaxis()->SetNdivisions(505);
            // TICKS X Axis
            geStat[2-j][i]->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
            geStat[2-j][i]->SetTitle("");
            geStat[2-j][i]->Draw("AP");
            gUrQMD[2-j][i]->Draw("l3");
            geUrQMDMirror[2-j][1-i]->Draw("l3");
            //geSys[icent][ew]->SetFillColorAlpha(1,0.5);
            //geSys[icent][ew]->Draw("3");
            geSys[2-j][i]->Draw("PS;;5 s=0.5");
            geSysMirror[2-j][1-i]->Draw("PS;;5 s=0.5");
            geStatMirror[2-j][1-i]->Draw("P");
            geStat[2-j][i]->Draw("P");
            l1[i]->Draw("same");
            l2[i]->Draw("same");
        }
    }
    C->cd();
    pad[0][2]->cd();
    le->Draw("same");
    C->SaveAs(OutputPDFName); 
    
    TCanvas *C1 = (TCanvas*) gROOT->FindObject("C1");
    if (C1) delete C1;
    C1 = new TCanvas("C1","canvas",800,1200);
    C1->SetFillStyle(4000);
    // Number of PADS
    const Int_t Nx1 = 2;
    const Int_t Ny1 = 3;
    // Canvas setup
    CanvasPartition(C1,Nx1,Ny1,lMargin,rMargin,bMargin,tMargin);
    TPad *pad1[Nx1][Ny1];
    for (Int_t i=0;i<Nx1;i++) {
        for (Int_t j=0;j<Ny1;j++) {
            C1->cd(0);
            // Get the pads previously created.
            char pname[16];
            sprintf(pname,"pad_%i_%i",i,j);
            pad1[i][j] = (TPad*) gROOT->FindObject(pname);
            pad1[i][j]->Draw();
            pad1[i][j]->SetFillStyle(4000);
            pad1[i][j]->SetFrameFillStyle(4000);
            pad1[i][j]->cd();
            // Size factors
            Float_t xFactor = pad1[0][0]->GetAbsWNDC()/pad1[i][j]->GetAbsWNDC();
            Float_t yFactor = pad1[0][0]->GetAbsHNDC()/pad1[i][j]->GetAbsHNDC();
            // Format for y axis
            geStat[2-j][i]->GetYaxis()->SetLabelFont(43);
            geStat[2-j][i]->GetYaxis()->SetLabelSize(24);
            geStat[2-j][i]->GetYaxis()->SetLabelOffset(0.02);
            geStat[2-j][i]->GetYaxis()->SetTitleFont(43);
            geStat[2-j][i]->GetYaxis()->SetTitleSize(24);
            geStat[2-j][i]->GetYaxis()->SetTitleOffset(2.5);
            //geStat[2-j][i]->GetYaxis()->CenterTitle();
            geStat[2-j][i]->GetYaxis()->SetNdivisions(505);
            // TICKS Y Axis
            geStat[2-j][i]->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
            // Format for x axis
            geStat[2-j][i]->GetXaxis()->SetLabelFont(43);
            geStat[2-j][i]->GetXaxis()->SetLabelSize(24);
            geStat[2-j][i]->GetXaxis()->SetLabelOffset(0.02);
            geStat[2-j][i]->GetXaxis()->SetTitleFont(43);
            geStat[2-j][i]->GetXaxis()->SetTitleSize(24);
            geStat[2-j][i]->GetXaxis()->SetTitleOffset(2.5);
            geStat[2-j][i]->GetXaxis()->CenterTitle();
            geStat[2-j][i]->GetXaxis()->SetNdivisions(505);
            // TICKS X Axis
            geStat[2-j][i]->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
            geStat[2-j][i]->SetTitle("");
            geStat[2-j][i]->Draw("AP");
            gUrQMD[2-j][i]->Draw("l3");
            //geSys[icent][ew]->SetFillColorAlpha(1,0.5);
            //geSys[icent][ew]->Draw("3");
            geSys[2-j][i]->Draw("PS;;5 s=0.5");
            //geSysMirror[2-j][1-i]->Draw("PS;;5 s=0.5");
            //geStatMirror[2-j][1-i]->Draw("P");
            geStat[2-j][i]->Draw("P");
            l1[i]->Draw("same");
            l2[i]->Draw("same");
        }
    }
    C1->cd();
    pad1[0][2]->cd();
    le->Draw("same");
    C1->SaveAs(OutputPDFName); 
    */
    // Margins
    Float_t lMargin3 = 0.2;
    Float_t rMargin3 = 0.05;
    Float_t bMargin3 = 0.15;
    Float_t tMargin3 = 0.05;
    /*
    TCanvas *C3 = (TCanvas*) gROOT->FindObject("C3");
    if (C3) delete C3;
    C3 = new TCanvas("C3","canvas",450,1200);
    C3->SetFillStyle(4000);
    // Number of PADS
    const Int_t Nx3 = 1;
    const Int_t Ny3 = 3;
    // Canvas setup
    CanvasPartition(C3,Nx3,Ny3,lMargin3,rMargin3,bMargin3,tMargin3);
    TPad *pad3[Nx3][Ny3];
    for (Int_t i=0;i<Nx3;i++) {
        for (Int_t j=0;j<Ny3;j++) {
            C3->cd(0);
            // Get the pads previously created.
            char pname[16];
            sprintf(pname,"pad_%i_%i",i,j);
            pad3[i][j] = (TPad*) gROOT->FindObject(pname);
            pad3[i][j]->Draw();
            pad3[i][j]->SetFillStyle(4000);
            pad3[i][j]->SetFrameFillStyle(4000);
            pad3[i][j]->cd();
            // Size factors
            Float_t xFactor = pad3[0][0]->GetAbsWNDC()/pad3[i][j]->GetAbsWNDC();
            Float_t yFactor = pad3[0][0]->GetAbsHNDC()/pad3[i][j]->GetAbsHNDC();
            // Format for y axis
            geSys[2-j][0]->GetYaxis()->SetLabelFont(43);
            geSys[2-j][0]->GetYaxis()->SetLabelSize(24);
            geSys[2-j][0]->GetYaxis()->SetLabelOffset(0.02);
            geSys[2-j][0]->GetYaxis()->SetTitleFont(43);
            geSys[2-j][0]->GetYaxis()->SetTitleSize(24);
            geSys[2-j][0]->GetYaxis()->SetTitleOffset(2);
            geSys[2-j][0]->GetYaxis()->CenterTitle();
            geSys[2-j][0]->GetYaxis()->SetNdivisions(505);
            // TICKS Y Axis
            geSys[2-j][0]->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
            // Format for x axis
            geSys[2-j][0]->GetXaxis()->SetLabelFont(43);
            geSys[2-j][0]->GetXaxis()->SetLabelSize(24);
            geSys[2-j][0]->GetXaxis()->SetLabelOffset(0.02);
            geSys[2-j][0]->GetXaxis()->SetTitleFont(43);
            geSys[2-j][0]->GetXaxis()->SetTitleSize(24);
            geSys[2-j][0]->GetXaxis()->SetTitleOffset(2.5);
            geSys[2-j][0]->GetXaxis()->CenterTitle();
            geSys[2-j][0]->GetXaxis()->SetNdivisions(505);
            // TICKS X Axis
            geSys[2-j][0]->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
            geSys[2-j][0]->GetXaxis()->SetLimits(-4.9,4.9); 
            geSys[2-j][0]->SetTitle(" "); 
            geSys[2-j][0]->Draw("APS;;5 s=0.5");
            gUrQMD[2-j][0]->Draw("l3");
            gUrQMD[2-j][1]->Draw("l3");
            gcos[2-j]->Draw("l3");
            //geSys[icent][ew]->SetFillColorAlpha(1,0.5);
            //geSys[icent][ew]->Draw("3");
            geSys[2-j][0]->Draw("PS;;5 s=0.5");
            geSys[2-j][1]->Draw("PS;;5 s=0.5");
            geStat[2-j][0]->Draw("P");
            geStat[2-j][1]->Draw("P");
            l2[0]->Draw("same");
            l2[1]->Draw("same");
            l3->Draw("same");
            l4->Draw("same");
        }
    }
    C3->cd();
    pad3[0][1]->cd();
    le2->Draw("same");
    C3->SaveAs(OutputPDFName); 
    */
    TCanvas *C4 = (TCanvas*) gROOT->FindObject("C4");
    if (C4) delete C4;
    C4 = new TCanvas("C4","canvas",450,1200);
    C4->SetFillStyle(4000);
    // Number of PADS
    const Int_t Nx4 = 1;
    //const Int_t Ny4 = 4;
    const Int_t Ny4 = 4;
    // Canvas setup
    CanvasPartition(C4,Nx4,Ny4,lMargin3,rMargin3,bMargin3,tMargin3);
    TPad *pad4[Nx4][Ny4];
    for (Int_t i=0;i<Nx4;i++) {
        for (Int_t j=0;j<Ny4;j++) {
            C4->cd(0);
            // Get the pads previously created.
            char pname[16];
            sprintf(pname,"pad_%i_%i",i,j);
            pad4[i][j] = (TPad*) gROOT->FindObject(pname);
            pad4[i][j]->Draw();
            pad4[i][j]->SetFillStyle(4000);
            pad4[i][j]->SetFrameFillStyle(4000);
            pad4[i][j]->cd();
            // Size factors
            Float_t xFactor = pad4[0][0]->GetAbsWNDC()/pad4[i][j]->GetAbsWNDC();
            Float_t yFactor = pad4[0][0]->GetAbsHNDC()/pad4[i][j]->GetAbsHNDC();
            // Format for y axis
            geStat[Ny4-1-j][0]->GetYaxis()->SetLabelFont(43);
            geStat[Ny4-1-j][0]->GetYaxis()->SetLabelSize(24);
            geStat[Ny4-1-j][0]->GetYaxis()->SetLabelOffset(0.02);
            geStat[Ny4-1-j][0]->GetYaxis()->SetTitleFont(43);
            geStat[Ny4-1-j][0]->GetYaxis()->SetTitleSize(24);
            geStat[Ny4-1-j][0]->GetYaxis()->SetTitleOffset(2.5);
            //geStat[Ny4-1-j][0]->GetYaxis()->CenterTitle();
            geStat[Ny4-1-j][0]->GetYaxis()->SetNdivisions(505);
            // TICKS Y Axis
            geStat[Ny4-1-j][0]->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
            // Format for x axis
            geStat[Ny4-1-j][0]->GetXaxis()->SetLabelFont(43);
            geStat[Ny4-1-j][0]->GetXaxis()->SetLabelSize(24);
            geStat[Ny4-1-j][0]->GetXaxis()->SetLabelOffset(0.02);
            geStat[Ny4-1-j][0]->GetXaxis()->SetTitleFont(43);
            geStat[Ny4-1-j][0]->GetXaxis()->SetTitleSize(24);
            geStat[Ny4-1-j][0]->GetXaxis()->SetTitleOffset(2.5);
            geStat[Ny4-1-j][0]->GetXaxis()->CenterTitle();
            geStat[Ny4-1-j][0]->GetXaxis()->SetNdivisions(505);
            // TICKS X Axis
            geStat[Ny4-1-j][0]->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
            geStat[Ny4-1-j][0]->SetTitle(" "); 
            geStat[Ny4-1-j][i]->Draw("AP");
            //gUrQMD[Ny4-1-j][i]->Draw("l3");
            //geUrQMDMirror[Ny4-1-j][1-i]->Draw("l3");
            //geSys[Ny4-1-j][i]->Draw("PS;;5 s=0.5");
            //geSysMirror[Ny4-1-j][1-i]->Draw("PS;;5 s=0.5");
            geStatMirror[Ny4-1-j][1-i]->Draw("P");
            geStat[Ny4-1-j][i]->Draw("P");
            l1[i]->Draw("same");
            l2[i]->Draw("same");
        }
    }
    C4->cd();
    pad4[0][1]->cd();
    le->Draw("same");
    C4->SaveAs(OutputPDFName);
    /*
    TCanvas *C2 = (TCanvas*) gROOT->FindObject("C2");
    if (C2) delete C2;
    C2 = new TCanvas("C2","canvas",1300,400);
    C2->SetFillStyle(4000);
    // Number of PADS
    const Int_t Nx2 = 3;
    const Int_t Ny2 = 1;
    // Canvas setup
    CanvasPartition(C2,Nx2,Ny2,lMargin,rMargin,bMargin,tMargin);
    TPad *pad2[Nx2][Ny2];
    for (Int_t i=0;i<Nx2;i++) {
        for (Int_t j=0;j<Ny2;j++) {
            C2->cd(0);
            // Get the pads previously created.
            char pname[16];
            sprintf(pname,"pad_%i_%i",i,j);
            pad2[i][j] = (TPad*) gROOT->FindObject(pname);
            pad2[i][j]->Draw();
            pad2[i][j]->SetFillStyle(4000);
            pad2[i][j]->SetFrameFillStyle(4000);
            pad2[i][j]->cd();
            
            // Size factors
            Float_t xFactor = pad2[0][0]->GetAbsWNDC()/pad2[i][j]->GetAbsWNDC();
            Float_t yFactor = pad2[0][0]->GetAbsHNDC()/pad2[i][j]->GetAbsHNDC();
            // Format for y axis
            geSys[i][0]->GetYaxis()->SetLabelFont(43);
            geSys[i][0]->GetYaxis()->SetLabelSize(24);
            geSys[i][0]->GetYaxis()->SetLabelOffset(0.02);
            geSys[i][0]->GetYaxis()->SetTitleFont(43);
            geSys[i][0]->GetYaxis()->SetTitleSize(24);
            geSys[i][0]->GetYaxis()->SetTitleOffset(1);
            geSys[i][0]->GetYaxis()->CenterTitle();
            geSys[i][0]->GetYaxis()->SetNdivisions(505);
            // TICKS Y Axis
            geSys[i][0]->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
            // Format for x axis
            geSys[i][0]->GetXaxis()->SetLabelFont(43);
            geSys[i][0]->GetXaxis()->SetLabelSize(24);
            geSys[i][0]->GetXaxis()->SetLabelOffset(0.02);
            geSys[i][0]->GetXaxis()->SetTitleFont(43);
            geSys[i][0]->GetXaxis()->SetTitleSize(24);
            geSys[i][0]->GetXaxis()->SetTitleOffset(1);
            geSys[i][0]->GetXaxis()->CenterTitle();
            geSys[i][0]->GetXaxis()->SetNdivisions(505);
            // TICKS X Axis
            geSys[i][0]->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
            geSys[i][0]->GetXaxis()->SetLimits(-4.9,4.9); 
            
            //geStat[i][0]->SetTitle("");
            geSys[i][0]->Draw("APS;;5 s=0.5");
            gcos[i]->Draw("l3");
            gUrQMD[i][0]->Draw("l3");
            gUrQMD[i][1]->Draw("l3");
            //geSys[icent][ew]->SetFillColorAlpha(1,0.5);
            //geSys[icent][ew]->Draw("3");
            geSys[i][0]->Draw("PS;;5 s=0.5");
            geSys[i][1]->Draw("PS;;5 s=0.5");
            geStat[i][0]->Draw("P");
            geStat[i][1]->Draw("P");
            l2[0]->Draw("same");
            l2[1]->Draw("same");
            l3->Draw("same");
            l4->Draw("same");
        }
    }
    C2->cd();
    pad2[1][0]->cd();
    le2->Draw("same");
    C2->SaveAs(OutputPDFName); 
    */
    /*
    //Compare with PHOBOS
    
    //v1 from PHOBOS vs. eta-ybeam  
    double x1[9]={1.84177,1.09264,0.34614,-0.358277,-0.806651,-1.30043,-1.79599,-2.29299,-2.79285};//19.6GeV
    double xerr1[9]={0.131344,0.132469,0.114251,0.0899632,0.101713,0.102112,0.101824,0.10216,0.102093};
    double y1[9]={0.119271,0.0912352,0.031701,-0.0139525,-0.0237081,-0.0191429,-0.013221,-0.00971628,-0.00185021};
    double yerr1[9]={0.0237411,0.0118685,0.00651729,0.00275923,0.00209791,0.0019658,0.00185551,0.00185868,0.00192151};

    double x2[9]={0.696226,-0.0601835,-0.806158,-1.50516,-1.95811,-2.45455,-2.95164,-3.45127,-3.95248};//62.4GeV
    double xerr2[9]={0.139005,0.130976,0.114685,0.0916473,0.102079,0.102118,0.10185,0.102171,0.102173};
    double y2[9]={0.0515044,0.00392858,-0.0169407,-0.0196981,-0.015187,-0.0124514,-0.00742903,-0.00345183,-0.00204235};
    double yerr2[9]={0.00571857,0.0011858,0.00119817,0.00164245,0.00127122,0.00133189,0.00105765,0.000955414,0.00106761};

    double x3[9]={-0.0311022,-0.787324,-1.53302,-2.24111,-2.69238,-3.18993,-3.6874,-4.1874,-4.68929};//130GeV
    double xerr3[9]={0.137215,0.13196,0.11519,0.0916862,0.102149,0.102202,0.101852,0.102119,0.102238};
    double y3[9]={0.00538043,-0.0127847,-0.0154049,-0.0209298,-0.0111316,-0.00764066,-0.000349314,-0.00114354,0.00263199};
    double yerr3[9]={0.00189709,0.00142972,0.00128735,0.00221135,0.0017566,0.00171018,0.00169515,0.00177189,0.00187604};

    double x4[9]={-0.452968,-1.21111,-1.95918,-2.66775,-3.11919,-3.61704,-4.11457,-4.61465,-5.11708};//200GeV
    double xerr4[9]={0.1394,0.131229,0.115226,0.0919191,0.102142,0.102117,0.101832,0.102162,0.102221};
    double y4[9]={-0.00106458,-0.007684,-0.0104365,-0.013613,-0.00802953,-0.0030862,-0.00244326,-0.0021339,0.000402647};
    double yerr4[9]={0.00160526,0.00123547,0.00104424,0.00152295,0.00119697,0.00132645,0.00129838,0.00145424,0.00183304};

    double x5[20]={-4.35,-4.25,-4.15,-4.05,-3.95,-3.85,-3.75,-3.65,-3.55,-3.45,-3.35,-3.25,-3.15,-3.05,-2.95,-2.85,-2.75,-2.65,-2.55,-2.45};
    double y5[20]={0.00989017,0.008895613,0.00787426,0.00675938,0.005880935,0.0046978,0.003549783,0.002552194,0.001584345,0.000632488,-0.00044001,-0.001598353,-0.002784383,-0.003921768,-0.004739398,-0.0056548,-0.006781605,-0.007762945,-0.00847943,-0.009728823};
    double xerr5[20]={0.0};
    double yerr5[20]={0.000191545,0.000184473,0.000181947,0.000181317,0.000181024,0.000180742,0.000180766,0.000180993,0.00018125,0.000181431,0.000181606,0.000181798,0.000181968,0.000182095,0.000182123,0.000182397,0.000182654,0.000183339,0.000186524,0.000195218};

    TGraphErrors* gv1ybeamStat = (TGraphErrors*)fStat->Get(Form("v1_vs_eta_minus_ybeam_charged_040"));
    TGraphErrors* gv1ybeamSys = (TGraphErrors*)fSys->Get(Form("v1_vs_eta_minus_ybeam_ch1_040_sys_error"));

    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.15);

    double xbg[2]={-6.0,3.0};
    double ybg[2]={-0.05,0.16};
    TGraph* bkgd = new TGraph(2,xbg,ybg);
    bkgd->Draw("AP");

    TH1D* tbg;
    TGraphErrors* gr1;
    TGraphErrors* gr2;
    TGraphErrors* gr3;
    TGraphErrors* gr4;
    TGraphErrors* gr5;
    
    tbg = new TH1D("t1","t1",35, -5.6,2.2);
    tbg->GetYaxis()->SetRangeUser(-0.05,0.18);
    tbg->SetStats(kFALSE);
    tbg->SetTitle("");
    tbg->GetYaxis()->SetTitle("v_{1}");
    tbg->GetXaxis()->SetTitle("#eta-y_{beam}");
    tbg->SetLineColor(1);
    tbg->GetYaxis()->SetLabelFont(43);
    tbg->GetYaxis()->SetLabelSize(24);
    tbg->GetYaxis()->SetLabelOffset(0.02);
    tbg->GetYaxis()->SetTitleFont(43);
    tbg->GetYaxis()->SetTitleSize(24);
    tbg->GetYaxis()->SetTitleOffset(1.5);
    tbg->GetYaxis()->CenterTitle();
    tbg->GetXaxis()->SetLabelFont(43);
    tbg->GetXaxis()->SetLabelSize(24);
    tbg->GetXaxis()->SetLabelOffset(0.02);
    tbg->GetXaxis()->SetTitleFont(43);
    tbg->GetXaxis()->SetTitleSize(24);
    tbg->GetXaxis()->SetTitleOffset(1.5);
    tbg->GetXaxis()->CenterTitle();
    tbg->Draw("");

    gr1 = new TGraphErrors(9,x1,y1,xerr1,yerr1);
    gr1->SetMarkerStyle(23);
    gr1->SetMarkerSize(1.5);
    gr1->SetMarkerColorAlpha(891,1);
    gr1->SetLineColorAlpha(891,1);
    gr1->Draw("same,P");
    
    gr2 = new TGraphErrors(9,x2,y2,xerr2,yerr2);
    gr2->SetMarkerStyle(22);
    gr2->SetMarkerSize(1.5);
    gr2->SetMarkerColorAlpha(601,1);
    gr2->SetLineColorAlpha(601,1);
    gr2->Draw("same,P");

    gr3 = new TGraphErrors(9,x3,y3,xerr3,yerr3);
    gr3->SetMarkerStyle(21);
    gr3->SetMarkerSize(1.5);
    gr3->SetMarkerColorAlpha(419,1);
    gr3->SetLineColorAlpha(419,1);
    gr3->Draw("same,P");
    
    gr4 = new TGraphErrors(9,x4,y4,xerr4,yerr4);
    gr4->SetMarkerStyle(33);
    gr4->SetMarkerSize(1.5);
    gr4->SetMarkerColorAlpha(882,1);
    gr4->SetLineColorAlpha(882,1);
    gr4->Draw("same,P");
    
    gr5 = new TGraphErrors(20,x5,y5,xerr5,yerr5);
    gr5->SetMarkerStyle(34);
    gr5->SetMarkerSize(1.5);
    gr5->SetMarkerColor(1);
    gr5->SetLineColor(1);
    gr5->Draw("same,P");

    for(int ieta=0;ieta<16;ieta++){
        gv1ybeamSys->SetPointError(ieta,0.15,gv1ybeamSys->GetErrorY(ieta));
    }

    gv1ybeamSys->SetMarkerStyle(20);
    gv1ybeamSys->SetMarkerSize(1.5);
    gv1ybeamSys->SetMarkerColorAlpha(1,0.8);
    gv1ybeamSys->SetLineColorAlpha(1,0.8);
    gv1ybeamSys->Draw("PS;;5 s=0.5");
 
    gv1ybeamStat->SetMarkerStyle(20);
    gv1ybeamStat->SetMarkerSize(1.5);
    gv1ybeamStat->SetMarkerColorAlpha(1,0.8);
    gv1ybeamStat->SetLineColorAlpha(1,0.8);
    gv1ybeamStat->Draw("same,P");

    TLegend* le3 = new TLegend(0.21,0.58,0.4,0.89);
    le3->AddEntry(gr1,"PHOBOS 19GeV","P");
    le3->AddEntry(gr2,"PHOBOS 62GeV","P");
    le3->AddEntry(gr3,"PHOBOS 130GeV","P");
    le3->AddEntry(gr4,"PHOBOS 200GeV","P");
    le3->AddEntry(gr5,"STAR 27GeV BESI","P");
    le3->AddEntry(gv1ybeamStat,"STAR 27GeV BESII","P");
    le3->SetBorderSize(0);
    le3->SetTextSize(0.03);
    le3->Draw("same");
    
    c1->SaveAs(OutputPDFName); 
    */
    TitlePage->SaveAs(Form("%s%s",OutputPDFName.Data(),"]"));//cannot .reset() for now because it will be saved again at the end
}
