void DrawSysCheck_ver2(TString mOutputFileName="DrawSysCheck_RCFrun27_tbt_ver2"){
    double labSize(0.06),titSize(0.08);
    gStyle->SetPalette(kBlueGreenYellow);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetTitleAlign(23);  // first number is horizontal 1=right/2=center/3=left  second number is vertical
    gStyle->SetLabelSize(0.05,"xy");
    gStyle->SetTitleSize(0.08,"xy");
    gStyle->SetTitleSize(0.08,"t");
    gStyle->SetTitleOffset(0.55,"xy");
    gROOT->ForceStyle();
    
    //TFile* fmeasuredv1 = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_code/v1NonFlow/v1_TGraphErrors_allitpc.root");
    TFile* fmeasuredv1 = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/root/tbt_v1_run26_combineVz_ver2.root");
    //Add measured v1 to this set of plots
    TGraphErrors* mv1measured[9];
    for(int icent=0;icent<9;icent++){
        //mv1measured[icent]=(TGraphErrors*)fmeasuredv1->Get(Form("v1Cent%dITPC%d",icent,0));
        mv1measured[icent]=(TGraphErrors*)fmeasuredv1->Get(Form("v1Cent%d",icent));
        mv1measured[icent]->SetMarkerStyle(21);
        mv1measured[icent]->SetMarkerSize(0.7);
        mv1measured[icent]->SetMarkerColorAlpha(882,1);
        mv1measured[icent]->SetLineColorAlpha(882,1);
        mv1measured[icent]->SetTitle("v_{1};#eta;v_{1}");
    }

    TFile* f[5];
    //default: mate dndeta, no pT dependent
    f[0] = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/root/DrawCorrectedv1_RCFrun26_tbt_default_Dec062022_ver2.root","READ");
    //fake phobos 27, no pT dependent 
    f[1] = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/root/DrawCorrectedv1_RCFrun26_tbt_dNdeta2_Dec062022_ver2.root","READ");
    //v1(pTsqrt)
    f[2] = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/root/DrawCorrectedv1_RCFrun26_tbt_pTsqrt_Dec062022_ver2.root","READ");
    //v1(pTsq)
    f[3] = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/root/DrawCorrectedv1_RCFrun26_tbt_pTsq_Dec062022_ver2.root","READ");
    //default rotate 180 degrees
    f[4] = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/root/DrawCorrectedv1_RCFrun26_tbt_Default_rot180_ver2.root","READ");

    TString OutputPDFName = "/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/pdf/"+mOutputFileName+".pdf";
    unique_ptr<TCanvas> TitlePage(new TCanvas("Title","Title",600,800));
    TPaveText tptitle(0.1,0.1,0.8,0.8);
    tptitle.AddText(Form("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/ver2/DrawSysCheck_ver2.C"));
    tptitle.Draw();
    TitlePage->SaveAs(Form("%s%s",OutputPDFName.Data(),"("));//cannot .reset() for now because it will be saved again at the end
   
    TGraphErrors* v1Corrected[5][9][2];//file, cent, all/charged particles 
    TGraphErrors* v1ybeam[5][9][2];//file, cent, all/charged particles

    //int color[7]={1,801,616,601,419,891,882};
    int color[7]={1,801,616,601,882};
    int centint[10]={80,70,60,50,40,30,20,10,5,0};

    for(int ifile=0;ifile<5;ifile++){
        for(int icent=2;icent<9;icent++){
            for(int ich=1;ich<2;ich++){
                v1Corrected[ifile][icent][ich]=(TGraphErrors*)f[ifile]->Get(Form("v1CorrectedCh%dCent%d",ich,icent));
                v1Corrected[ifile][icent][ich]->SetMarkerStyle(20);
                v1Corrected[ifile][icent][ich]->SetMarkerSize(0.5);
                v1Corrected[ifile][icent][ich]->SetMarkerColorAlpha(color[ifile],1);
                v1Corrected[ifile][icent][ich]->SetLineColorAlpha(color[ifile],1);
                v1Corrected[ifile][icent][ich]->SetMaximum(0.17);
                v1Corrected[ifile][icent][ich]->SetMinimum(-0.17);
                v1Corrected[ifile][icent][ich]->GetXaxis()->SetTitle("#eta");
                v1Corrected[ifile][icent][ich]->GetYaxis()->SetTitle("v_{1}");
                v1Corrected[ifile][icent][ich]->SetTitle(Form("%d~%d%s",centint[icent+1],centint[icent],"%"));
            }
            v1ybeam[ifile][icent][0]=(TGraphErrors*)f[ifile]->Get(Form("v1_vs_eta_minus_ybeam_all_040"));
            v1ybeam[ifile][icent][1]=(TGraphErrors*)f[ifile]->Get(Form("v1_vs_eta_minus_ybeam_charged_040"));
        }
    }
    
    double range[9]={0.17,0.17,0.17,0.17,0.17,0.17,0.17,0.12,0.07}; 
    TLine* l1=new TLine(-5.2,0.0,5.2,0.0);
    TLine* l2[9];
    TLine* l3[9];
    TLine* l4[9];
    for(int icent=0;icent<9;icent++){
        l2[icent]=new TLine(0.0,-range[icent],0.0,range[icent]);
        l3[icent]=new TLine(-3.4,-range[icent],-3.4,range[icent]);
        l4[icent]=new TLine(3.4,-range[icent],3.4,range[icent]);
        l2[icent]->SetLineStyle(2);
        l3[icent]->SetLineStyle(2);
        l4[icent]->SetLineStyle(2);
    }
    l1->SetLineStyle(2);

    TLegend* lcom=new TLegend(0.1,0.15,0.3,0.85); 
    lcom->AddEntry(v1Corrected[0][2][1],Form("Unfolded #frac{dN}{d#eta}, no p_{T} dep."));
    lcom->AddEntry(v1Corrected[1][2][1],Form("Fake PHOBOS 27 GeV #frac{dN}{d#eta}, no p_{T} dep."));
    lcom->AddEntry(v1Corrected[2][2][1],Form("Unfolded #frac{dN}{d#eta}, quadratic p_{T}"));
    lcom->AddEntry(v1Corrected[3][2][1],Form("Unfolded #frac{dN}{d#eta}, p_{T} square"));
    //lcom->AddEntry(v1Corrected[4][2][1],Form("|#eta|<0.8 Mate #frac{dN}{d#eta}, no pT dep."));
    //lcom->AddEntry(v1Corrected[5][2][1],Form("|#eta|<0.6 Mate #frac{dN}{d#eta}, no pT dep."));
    lcom->AddEntry(v1Corrected[4][2][1],Form("Rotate 180 degree"));
    //lcom->AddEntry(mv1measured[2],Form("before correction"));
    lcom->SetBorderSize(0);
    lcom->SetTextSize(0.06);
    
    TCanvas* ccom = new TCanvas("ccom","ccom",1200,600);
    ccom->Divide(4,2,0,0);
    for(int icent=2;icent<9;icent++){
        if(icent<5) ccom->cd(icent-1);
        else ccom->cd(icent);
        v1Corrected[4][icent][1]->GetXaxis()->SetRangeUser(-5.2,5.2); 
        v1Corrected[4][icent][1]->GetYaxis()->CenterTitle();
        v1Corrected[4][icent][1]->GetXaxis()->CenterTitle();
        v1Corrected[4][icent][1]->Draw("AP");
        v1Corrected[1][icent][1]->Draw("P");
        v1Corrected[2][icent][1]->Draw("P");
        v1Corrected[3][icent][1]->Draw("P");
        //v1Corrected[5][icent][1]->Draw("P");
        //v1Corrected[6][icent][1]->Draw("P");
        v1Corrected[0][icent][1]->Draw("P");
        //mv1measured[icent]->Draw("P");
        l1->Draw("same");
        l2[0]->Draw("same");
        l3[0]->Draw("same");
        l4[0]->Draw("same");
    }
    ccom->cd(4);
    lcom->Draw("same");

    ccom->SaveAs(OutputPDFName); 

    double xInconsis[9][2][32][5];//cent, all/charged particles, 32 v1 points, file
    double yInconsis[9][2][32][5];
    double xerrInconsis[9][2][32][5];
    double yerrInconsis[9][2][32][5];

    double xAll[9][2][32][5];//cent, all/charged particles, 32 v1 points, file
    double yAll[9][2][32][5];
    double xerrAll[9][2][32][5];
    double yerrAll[9][2][32][5];
    
    double xSys[9][2][32];
    double ySys[9][2][32];
    double xerrSys[9][2][32];
    double yerrSys[9][2][32];

    for(int icent=2;icent<9;icent++){
        for(int ich=1;ich<2;ich++){
            for(int i=0;i<32;i++){
                xSys[icent][ich][i]=v1Corrected[0][icent][ich]->GetX()[i];
                xerrSys[icent][ich][i]=0.0;
                ySys[icent][ich][i]=v1Corrected[0][icent][ich]->GetY()[i];
                double tmp_sys_err=0.0;

                for(int ifile=0;ifile<5;ifile++){
                    xAll[icent][ich][i][ifile]=ifile+1;
                    xerrAll[icent][ich][i][ifile]=0.0; 
                    yAll[icent][ich][i][ifile]=v1Corrected[ifile][icent][ich]->GetY()[i];
                    yerrAll[icent][ich][i][ifile]=v1Corrected[ifile][icent][ich]->GetErrorY(i);
                    xInconsis[icent][ich][i][ifile]=ifile+1;
                    xerrInconsis[icent][ich][i][ifile]=0.0; 
                    if(ifile!=0){
                        double tmp_diff=fabs(v1Corrected[0][icent][ich]->GetY()[i]-v1Corrected[ifile][icent][ich]->GetY()[i]);
                        double tmp_differr=sqrt(fabs(pow(v1Corrected[0][icent][ich]->GetErrorY(i),2)-pow(v1Corrected[ifile][icent][ich]->GetErrorY(i),2)));
                        if(tmp_diff>tmp_differr){
                            yInconsis[icent][ich][i][ifile]=v1Corrected[ifile][icent][ich]->GetY()[i];
                            yerrInconsis[icent][ich][i][ifile]=v1Corrected[ifile][icent][ich]->GetErrorY(i);
                            tmp_sys_err+=fabs(pow(tmp_diff,2)-pow(tmp_differr,2))/12.0;
                            //tmp_sys_err+=fabs(pow(tmp_diff,2)-pow(tmp_differr,2));
                        }
                        else{
                            yInconsis[icent][ich][i][ifile]=999;
                            yerrInconsis[icent][ich][i][ifile]=0.0;
                        }
                    }
                    else{
                        yInconsis[icent][ich][i][ifile]=999;
                        yerrInconsis[icent][ich][i][ifile]=0.0;
                    }
                }
                yerrSys[icent][ich][i]=sqrt(tmp_sys_err);
            }
        }
    }

    
    TGraphErrors* GInconsis[9][2][32];
    TGraphErrors* GAll[9][2][32];
    TGraphErrors* GSysError[9][2];

    for(int icent=2;icent<9;icent++){
        for(int ich=1;ich<2;ich++){
            GSysError[icent][ich]=new TGraphErrors(32,xSys[icent][ich],ySys[icent][ich],xerrSys[icent][ich],yerrSys[icent][ich]); 
            GSysError[icent][ich]->SetMarkerStyle(20);
            GSysError[icent][ich]->SetMarkerSize(0.6);
            GSysError[icent][ich]->SetMarkerColor(1);
            
            for(int i=0;i<32;i++){
                GAll[icent][ich][i]=new TGraphErrors(5,xAll[icent][ich][i],yAll[icent][ich][i],xerrAll[icent][ich][i],yerrAll[icent][ich][i]); 
                GAll[icent][ich][i]->SetMarkerStyle(20);
                GAll[icent][ich][i]->SetMarkerSize(0.6);
                GAll[icent][ich][i]->SetMarkerColor(1);
                GAll[icent][ich][i]->SetLineColor(1);
                GAll[icent][ich][i]->SetTitle(Form("%d~%d%s #eta=%.2f;sys. check;v1(#eta)",centint[icent+1],centint[icent],"%",v1Corrected[0][icent][ich]->GetX()[i]));
                GInconsis[icent][ich][i]=new TGraphErrors(5,xInconsis[icent][ich][i],yInconsis[icent][ich][i],xerrInconsis[icent][ich][i],yerrInconsis[icent][ich][i]); 
                GInconsis[icent][ich][i]->SetMarkerStyle(20);
                GInconsis[icent][ich][i]->SetMarkerSize(0.6);
                GInconsis[icent][ich][i]->SetMarkerColor(2);
                GInconsis[icent][ich][i]->SetLineColor(2);
            }
        }
    }

    TCanvas* c[9][2];
    for(int icent=2;icent<9;icent++){
        for(int iw=0;iw<2;iw++){
            c[icent][iw]=new TCanvas(Form("cCent%dEW%d",icent,iw),Form("cCent%dEW%d",icent,iw),800,1500);
            c[icent][iw]->Divide(2,8);
        }
        for(int i=0;i<16;i++){
            c[icent][0]->cd(i+1);
            GAll[icent][1][i]->Draw("AP");
            GInconsis[icent][1][i]->Draw("P");
            c[icent][1]->cd(i+1);
            GAll[icent][1][i+16]->Draw("AP");
            GInconsis[icent][1][i+16]->Draw("P");
        }
        for(int iw=0;iw<2;iw++){
            c[icent][iw]->SaveAs(OutputPDFName); 
        }
    }
    
    TCanvas* cSys = new TCanvas("cSys","cSys",1200,600);
    cSys->Divide(4,2);
    for(int icent=2;icent<9;icent++){
        cSys->cd(icent);
        v1Corrected[0][icent][1]->SetMaximum(range[icent]); 
        v1Corrected[0][icent][1]->SetMinimum(-range[icent]); 
        v1Corrected[0][icent][1]->GetYaxis()->CenterTitle();
        v1Corrected[0][icent][1]->Draw("AP");
        GSysError[icent][1]->Draw("[]"); 
        l1->Draw("same");
        l2[icent]->Draw("same");
        l3[icent]->Draw("same");
        l4[icent]->Draw("same");
    }
    cSys->SaveAs(OutputPDFName); 

    TitlePage->SaveAs(Form("%s%s",OutputPDFName.Data(),"]"));//cannot .reset() for now because it will be saved again at the end

    //Combine centrality to compare with PHOBOS 
    double x040[2][16];
    double y040[2][16];
    double xerr040[2][16];
    double yerr040[2][16];

    for(int ich=1;ich<2;ich++){
        for(int iring=0;iring<16;iring++){

            double tmp_sum=0.0;
            double tmp_sum_err=0.0;
            
            for(int icent=8;icent>3;icent--){
                if(icent>6){
                    tmp_sum+=(-ySys[icent][ich][iring]+ySys[icent][ich][31-iring])/16.0;
                    tmp_sum_err+=pow(yerrSys[icent][ich][iring]/16.0,2)+pow(yerrSys[icent][ich][31-iring]/16.0,2);
                }
                else{
                    tmp_sum+=(-ySys[icent][ich][iring]+ySys[icent][ich][31-iring])/8.0;
                    tmp_sum_err+=pow(yerrSys[icent][ich][iring]/8.0,2)+pow(yerrSys[icent][ich][31-iring]/8.0,2);

                }
            }
            y040[ich][iring]=tmp_sum;
            yerr040[ich][iring]=sqrt(tmp_sum_err);
            xerr040[ich][iring]=0.0;
            x040[ich][iring]=xSys[8][ich][31-iring]-3.4;
        }

    }

    TGraphErrors* GSysError040[2];
    
    for(int ich=1;ich<2;ich++){
        GSysError040[ich] = new TGraphErrors(16,x040[ich],y040[ich],xerr040[ich],yerr040[ich]);
    }
    
    
    TString OutputROOTName = "/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/root/"+mOutputFileName+".root";
    TFile *foutput = new TFile(OutputROOTName.Data(),"RECREATE");
    for(int ich=1;ich<2;ich++){//0_all,1_charged particles 
        for(int icent=2;icent<9;icent++){
            GSysError[icent][ich]->SetName(Form("v1CorrectedSysErrorCh%dCent%d",ich,icent));
            GSysError[icent][ich]->Write(); 
        }
        GSysError040[ich]->SetName(Form("v1_vs_eta_minus_ybeam_ch%d_040_sys_error",ich));
        GSysError040[ich]->Write();
    }
    foutput->Write();
    foutput->Close();

}