void DrawCorrectedv1_ver2(TString mOutputFileName="DrawCorrectedv1_RCFrun26_tbt_pTsqrt_Dec062022_ver2"){
    double labSize(0.06),titSize(0.08);
    gStyle->SetTitleSize(titSize,"xyT");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetPadLeftMargin(0.17);
    gStyle->SetTitleAlign(23);  // first number is horizontal 1=right/2=center/3=left  second number is vertical
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetLabelSize(0.08,"xy");
    gStyle->SetTitleSize(0.08,"xy");
    gStyle->SetTitleSize(0.08,"t");
    gStyle->SetTitleOffset(0.4,"xy");
    TFile* f = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/pTsqrt/Sep212022_itr1_v1pTsqrtonTree.root","READ");
    //TFile* finputv1 = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/dNdeta2/MakeNewInput_dndetaFakePhobos_fitasy_RCFrun26_dNdeta2_itr2.root","READ");
    TFile* fmeasuredv1 = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/root/tbt_v1_run26_combineVz_ver2.root");

    TString OutputPDFName = "/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/pdf/"+mOutputFileName+".pdf";
    unique_ptr<TCanvas> TitlePage(new TCanvas("Title","Title",600,800));
    TPaveText tptitle(0.1,0.1,0.8,0.8);
    tptitle.AddText(Form("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/ver2/DrawCorrectedv1_ver2.C"));
    tptitle.AddText(Form("inputv1:/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/pTsqrt/MakeNewInput_dndetaMate_fitasy_RCFrun26_default_itr3.root"));
    tptitle.AddText(Form("GEANTv1:/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/pTsqrt/Sep212022_itr1_v1pTsqrtonTree.root"));
    tptitle.AddText(Form("measured v1:/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/root/tbt_v1_run26_combineVz_ver2.root"));
    tptitle.Draw();
    TitlePage->SaveAs(Form("%s%s",OutputPDFName.Data(),"("));//cannot .reset() for now because it will be saved again at the end

    int color[9]={801,632,891,616,886,601,861,417,419};
    int vzbin[10]={-45, -35, -25, -15, -5, 5, 15, 25, 35, 45};
    int centint[10]={80,70,60,50,40,30,20,10,5,0};
    
    /* 
    TH1D* v1Input[9];
    for(int icent=2;icent<9;icent++){
        v1Input[icent] =(TH1D*)finputv1->Get(Form("v1Cent%dHisto600",icent));
        v1Input[icent]->SetTitle("v_{1};#eta;v_{1}");
        v1Input[icent]->SetTitle(Form("%d~%d%s",centint[icent+1],centint[icent],"%"));
        v1Input[icent]->SetLineColor(4);
        v1Input[icent]->SetMarkerColor(4);
    }
    */
    TProfile* v1MC[2][9][9][2];//all/charged particles,cent, vz, weighting schemes  
    TProfile* v1EPD[9][9][2];
    TProfile2D* mAveEta[9][2];//cent, weighting schemes
    
    for(int icent=2;icent<9;icent++){
        for(int iw=0;iw<2;iw++){
            mAveEta[icent][iw]=(TProfile2D*)f->Get(Form("AveEtaCent%dWt%d",icent,iw));
            for(int ivz=0;ivz<9;ivz++){
                v1MC[0][icent][ivz][iw]=(TProfile*)f->Get(Form("v1McTracks_Cent%diz%dWt%d",icent,ivz,iw));
                v1MC[0][icent][ivz][iw]->SetTitle(Form("%d~%d%s %d<vz<%d cm",centint[icent+1],centint[icent],"%",vzbin[ivz],vzbin[ivz+1]));
                v1MC[0][icent][ivz][iw]->SetLineColor(color[ivz]);
                v1MC[1][icent][ivz][iw]=(TProfile*)f->Get(Form("v1McTracksCh_Cent%diz%dWt%d",icent,ivz,iw));
                v1MC[1][icent][ivz][iw]->SetTitle(Form("%d~%d%s %d<vz<%d cm",centint[icent+1],centint[icent],"%",vzbin[ivz],vzbin[ivz+1]));
                v1MC[1][icent][ivz][iw]->SetLineColor(color[ivz]);
                v1EPD[icent][ivz][iw]=(TProfile*)f->Get(Form("v1EPDhits_Cent%diz%dWt%d",icent,ivz,iw));
            }
        }
    }

    //Average v1MC over all vz bins 
    TH1D* Avev1MC[2][9];//all/charged particles, cent
    for(int ich=0;ich<2;ich++){
        for(int icent=2;icent<9;icent++){
            Avev1MC[ich][icent]=new TH1D(Form("Avev1MCCh%dCent%d",ich,icent),Form("Avev1MCCh%dCent%d",ich,icent),120,-6.0,6.0);
            for(int i=0;i<v1MC[1][2][0][1]->GetNbinsX();i++){
                double tmp=0.0;
                double tmp_err=0.0;
                for(int ivz=0;ivz<9;ivz++){
                    tmp+=v1MC[ich][icent][ivz][1]->GetBinContent(i+1);
                    tmp_err+=pow(v1MC[ich][icent][ivz][1]->GetBinError(i+1),2);
                }
                tmp/=9.0;
                tmp_err=sqrt(tmp_err)/9.0;
                Avev1MC[ich][icent]->SetBinContent(i+1,tmp);
                Avev1MC[ich][icent]->SetBinError(i+1,tmp_err);
            }
        }
    } 

    //Add measured v1 to this set of plots
    TGraphErrors* mv1measured[9];
    for(int icent=0;icent<9;icent++){
        //mv1measured[icent]=(TGraphErrors*)fmeasuredv1->Get(Form("v1Cent%dITPC%d",icent,ITPC));
        mv1measured[icent]=(TGraphErrors*)fmeasuredv1->Get(Form("v1Cent%d",icent));
        mv1measured[icent]->SetMarkerStyle(21);
        mv1measured[icent]->SetMarkerSize(0.7);
        mv1measured[icent]->SetMarkerColorAlpha(882,1);
        mv1measured[icent]->SetLineColorAlpha(882,1);
        mv1measured[icent]->SetTitle("v_{1};#eta;v_{1}");
        mv1measured[icent]->SetTitle(Form("%d~%d%s",centint[icent+1],centint[icent],"%"));
    }
    TFile* furqmd = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_code/UrQMD/UrQMDAna/Feb122022_run11_UrQMD_v1Hists.root","READ");
    //Corresponding QA file for the UrQMD data
    TFile* fQA = new TFile("/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_code/UrQMD/UrQMDAna/Jan032022_run2_UrQMD_QAHists.root","READ");
    
    TProfile* v1EPRaw[9][2];
    TProfile* mResolution[9];
    TProfile* mCosphi[9];
    
    double binsetaW[17]={2.14,2.2,2.27,2.34,2.41,2.5,2.59,2.69,2.81,2.94,3.08,3.26,3.47,3.74,4.03,4.42,5.09};
    double binseta[2][17];
    for(int ibin=0;ibin<17;ibin++){
        binseta[0][ibin]=-binsetaW[ibin];
        binseta[1][ibin]=binsetaW[ibin];
    }

    for(int icent=0;icent<9;icent++){
        mResolution[icent]=(TProfile*)furqmd->Get(Form("ResolutionTPC%dCent%dFF%d",13,icent,0));
        for(int ew=0;ew<2;ew++){
            v1EPRaw[icent][ew]=(TProfile*)furqmd->Get(Form("v1EPDEPRawCent%dEW%dREF%dFF%d",icent,ew,13,0));
        }
        mCosphi[icent]=(TProfile*)fQA->Get(Form("CosphiCent%d",icent));
        mCosphi[icent]->SetLineColor(2);
        mCosphi[icent]->SetMarkerColor(2);
    }
    
    double xurqmdv1[9][2][16];//cent,ew,ring
    double yurqmdv1[9][2][16];//cent,ew,ring
    double xerrurqmdv1[9][2][16];//cent,ew,ring
    double yerrurqmdv1[9][2][16];//cent,ew,ring
    
    for(int cent=0;cent<9;cent++){
        for(int ew=0;ew<2;ew++){
            double tmp_res=sqrt(mResolution[cent]->GetBinContent(1)*mResolution[cent]->GetBinContent(2)/mResolution[cent]->GetBinContent(3));
            for(int ir=0;ir<16;ir++){
                double tmp_eta=0.5*(binseta[ew][ir]+binseta[ew][ir+1]);
                double ringetabin=v1EPRaw[cent][ew]->FindBin(tmp_eta);
            
                xurqmdv1[cent][ew][ir]=tmp_eta;
                yurqmdv1[cent][ew][ir]=v1EPRaw[cent][ew]->GetBinContent(ringetabin)/tmp_res;
                xerrurqmdv1[cent][ew][ir]=0.0;
                yerrurqmdv1[cent][ew][ir]=v1EPRaw[cent][ew]->GetBinError(ringetabin)/tmp_res;
            }
        }
    }
    TGraphErrors* Gurqmdv1[9][2];
    for(int cent=0;cent<9;cent++){
        for(int ew=0;ew<2;ew++){
            Gurqmdv1[cent][ew]=new TGraphErrors(16,xurqmdv1[cent][ew],yurqmdv1[cent][ew],xerrurqmdv1[cent][ew],yerrurqmdv1[cent][ew]);
            Gurqmdv1[cent][ew]->SetMarkerStyle(20);
            Gurqmdv1[cent][ew]->SetMarkerSize(0.8);
            Gurqmdv1[cent][ew]->SetMarkerColor(4);
            Gurqmdv1[cent][ew]->SetLineColor(4);
            Gurqmdv1[cent][ew]->GetXaxis()->SetTitle("#eta");
            Gurqmdv1[cent][ew]->GetYaxis()->SetTitle("v1(#eta)");
            //Gurqmdv1[cent][ew][itpc]->SetTitle(Form("%d~%d%s",centint[cent+1],centint[cent],"%"));
        }
    }

    double x[9][9][32];//cent, vz, rings on both sides
    double y[9][9][32];
    double xerr[9][9][32];
    double yerr[9][9][32];

    //conbine 9 points in a group 
    double xCom[9][32][9];
    double yCom[9][32][9];
    double xerrCom[9][32][9]={0.0};
    double yerrCom[9][32][9];//cent,ieta, group,  ring
    
    double xCloud[9][32];//cent,9 groups
    double yCloud[9][32];  
    double xerrCloud[9][32]={0.0};  
    double yerrCloud[9][32];  

    double xv1InNew[2][9][32];//all/charged particles
    double yv1InNew[2][9][32];
    double xerrv1InNew[2][9][32];
    double yerrv1InNew[2][9][32];
    
    //v1(TPC), with respect to EPDE+W, resolution calculated by sqrt(2cos(PsiE-PsiW))
    double x1TPC[9][10]={-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000,-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000,-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000,-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000,-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000,-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000,-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000,-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000,-0.900000,-0.700000,-0.500000,-0.300000,-0.100000,0.100000,0.300000,0.500000,0.700000,0.900000};
    double y1TPC[9][10]={0.015377,0.010634,0.006487,0.002240,-0.000186,-0.002736,-0.005473,-0.009323,-0.012413,-0.017503,0.017499,0.012069,0.008149,0.004067,0.000951,-0.002547,-0.006373,-0.009553,-0.013584,-0.018322,0.017144,0.012519,0.008166,0.004470,0.001157,-0.002216,-0.005694,-0.009531,-0.013233,-0.017416,0.015656,0.011882,0.008051,0.004741,0.001444,-0.001724,-0.005360,-0.008639,-0.012363,-0.015997,0.013847,0.010607,0.007200,0.004348,0.001308,-0.001656,-0.004583,-0.007580,-0.010777,-0.013976,0.011218,0.008779,0.006168,0.003524,0.001112,-0.001296,-0.003870,-0.006311,-0.008826,-0.011489,0.008546,0.006652,0.004673,0.002796,0.000908,-0.001089,-0.002974,-0.004819,-0.006705,-0.008703,0.006226,0.004925,0.003295,0.001982,0.000665,-0.000830,-0.002201,-0.003488,-0.005074,-0.006444,0.004604,0.003681,0.002522,0.001268,0.000458,-0.000509,-0.001574,-0.002219,-0.003455,-0.004670};
    double xerr1TPC[9][10]={0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};
    double yerr1TPC[9][10]={0.000509,0.000497,0.000502,0.000510,0.000509,0.000509,0.000507,0.000499,0.000493,0.000499,0.000256,0.000250,0.000250,0.000252,0.000251,0.000251,0.000250,0.000249,0.000248,0.000251,0.000148,0.000144,0.000144,0.000145,0.000145,0.000144,0.000144,0.000143,0.000143,0.000145,0.000097,0.000095,0.000095,0.000095,0.000095,0.000095,0.000094,0.000094,0.000094,0.000095,0.000069,0.000067,0.000067,0.000067,0.000067,0.000067,0.000067,0.000067,0.000067,0.000068,0.000056,0.000055,0.000054,0.000055,0.000054,0.000054,0.000054,0.000054,0.000054,0.000055,0.000052,0.000051,0.000050,0.000050,0.000050,0.000050,0.000050,0.000050,0.000050,0.000051,0.000091,0.000088,0.000088,0.000088,0.000088,0.000087,0.000087,0.000087,0.000087,0.000089,0.000141,0.000135,0.000134,0.000134,0.000134,0.000134,0.000133,0.000133,0.000134,0.000137};

    for(int icent=2;icent<9;icent++){
        vector<tuple<double, double, double,double>> v;
        for(int ivz=0;ivz<9;ivz++){
            for(int i=0;i<32;i++){
                x[icent][ivz][i]=mAveEta[icent][1]->GetBinContent(i+1,ivz+1);
                y[icent][ivz][i]=v1EPD[icent][ivz][1]->GetBinContent(i+1);
                xerr[icent][ivz][i]=0.0;
                yerr[icent][ivz][i]=v1EPD[icent][ivz][1]->GetBinError(i+1);
                v.push_back(make_tuple(x[icent][ivz][i],y[icent][ivz][i],xerr[icent][ivz][i],yerr[icent][ivz][i]));
            }
        }
        sort(v.begin(),v.end());
        for(int ig=0;ig<32;ig++){
            for(int id=0;id<9;id++){//group every 9 points
                xCom[icent][ig][id]=get<0>(v[ig*9+id]);
                yCom[icent][ig][id]=get<1>(v[ig*9+id]);
                xerrCom[icent][ig][id]=get<2>(v[ig*9+id]);
                yerrCom[icent][ig][id]=get<3>(v[ig*9+id]);
            }
            double xsum=0.0;
            double ysum=0.0;
            double err2sum=0.0;
            for(int id=0;id<9;id++){
            xsum+=xCom[icent][ig][id];
            ysum+=yCom[icent][ig][id]; 
            err2sum+=pow(yerrCom[icent][ig][id],2);
            }
            xsum/=9.0;
            ysum/=9.0;
            err2sum=sqrt(err2sum)/9.0;
            xCloud[icent][ig]=xsum;
            yCloud[icent][ig]=ysum;
            yerrCloud[icent][ig]=err2sum;
            //cout<<yerrCloud[icent][ig]<<endl;

            /*
            xv1InNew[icent][ig]=xCloud[icent][ig];
            yv1InNew[icent][ig]=v1Input[icent]->GetBinContent(v1Input[icent]->FindBin(xCloud[icent][ig]))
            +mv1measured[icent]->GetY()[ig]-yCloud[icent][ig];
            xerrv1InNew[icent][ig]=0.0;
            //yerrv1InNew[icent][ig]=0.0;
            yerrv1InNew[icent][ig]=sqrt(pow(mv1measured[icent]->GetErrorY(ig),2)+pow(yerrCloud[icent][ig],2));
            //yerrv1InNew[icent][ig]=sqrt(pow(mv1measured[icent]->GetErrorY(ig),2));
            */
            for(int ich=0;ich<2;ich++){
                xv1InNew[ich][icent][ig]=xCloud[icent][ig];
                yv1InNew[ich][icent][ig]=Avev1MC[ich][icent]->GetBinContent(Avev1MC[ich][icent]->FindBin(xCloud[icent][ig]))
                +mv1measured[icent]->GetY()[ig]-yCloud[icent][ig];
                xerrv1InNew[ich][icent][ig]=0.0;
                yerrv1InNew[ich][icent][ig]=sqrt(pow(mv1measured[icent]->GetErrorY(ig),2)+pow(yerrCloud[icent][ig],2)+
                pow(Avev1MC[ich][icent]->GetBinError(Avev1MC[ich][icent]->FindBin(xCloud[icent][ig])),2));
            }
        }
        v.clear();
    }

    TGraphErrors* gv1Corrected[2][9];//corrected by v1(all particles), v1(charged particles)
    for(int ich=0;ich<2;ich++){
        for(int icent=2;icent<9;icent++){
            gv1Corrected[ich][icent]=new TGraphErrors(32,xv1InNew[ich][icent],yv1InNew[ich][icent],xerrv1InNew[ich][icent],yerrv1InNew[ich][icent]);
            gv1Corrected[ich][icent]->SetMarkerSize(0.7);
            gv1Corrected[ich][icent]->SetMarkerStyle(21);
            gv1Corrected[ich][icent]->SetMarkerColorAlpha(1,1);
            gv1Corrected[ich][icent]->SetLineColorAlpha(1,1);
            gv1Corrected[ich][icent]->GetXaxis()->SetTitle("#eta");
            gv1Corrected[ich][icent]->GetYaxis()->SetTitle("v1(#eta)");
            gv1Corrected[ich][icent]->SetTitle(Form("%d~%d%s",centint[icent+1],centint[icent],"%"));
        }
    }

      
    TLegend* le =new TLegend(0.3,0.05,0.5,0.35);
    le->AddEntry(mv1measured[2],"before correction","PL"); 
    le->AddEntry(gv1Corrected[0][2],"after correction,all","PL"); 
    le->AddEntry(gv1Corrected[1][2],"after correction,charged","PL"); 
    le->SetBorderSize(0);
    le->SetTextSize(0.04);
    
    TLegend* leurqmd=new TLegend(0.3,0.15,0.4,0.45); 
    leurqmd->AddEntry(gv1Corrected[0][2],Form("Data #Psi^{|#eta<1.0|},all"));
    leurqmd->AddEntry(gv1Corrected[1][2],Form("Data #Psi^{|#eta<1.0|},charged"));
    leurqmd->AddEntry(Gurqmdv1[2][0],Form("UrQMD #Psi^{|#eta<1.0|}"));
    leurqmd->AddEntry(mCosphi[2],"UrQMD RP=0","PL");
    leurqmd->SetBorderSize(0);
    leurqmd->SetTextSize(0.06);

    TLine* l1=new TLine(-5.2,0.0,5.2,0.0);
    TLine* l2=new TLine(0.0,-0.18,0.0,0.18);
    TLine* l3=new TLine(-3.4,-0.18,-3.4,0.18);
    TLine* l4=new TLine(3.4,-0.18,3.4,0.18);
    l1->SetLineStyle(2);
    l2->SetLineStyle(2);
    l3->SetLineStyle(2);
    l4->SetLineStyle(2);

    double max_array[9]={0.0,0.0,0.28,0.26,0.26,0.22,0.18,0.11,0.06};
    TCanvas* c[9]; 
    for(int icent=2;icent<9;icent++){
        c[icent]=new TCanvas(Form("c%d",icent),Form("c%d",icent),800,600);
        c[icent]->cd();
        gv1Corrected[1][icent]->SetMarkerColor(3);
        gv1Corrected[1][icent]->SetLineColor(3);
        gv1Corrected[0][icent]->Draw("AP");
        gv1Corrected[1][icent]->Draw("P");
        mv1measured[icent]->Draw("P");
        le->Draw("same");
        c[icent]->SaveAs(OutputPDFName);
    }
    
    TLegend* lcom=new TLegend(0.3,0.15,0.4,0.45); 
    lcom->AddEntry(gv1Corrected[0][2],Form("After correction, all"));
    lcom->AddEntry(gv1Corrected[1][2],Form("After correction, charged"));
    lcom->AddEntry(mv1measured[2],Form("Before correction"));
    lcom->SetBorderSize(0);
    lcom->SetTextSize(0.06);
    
    TCanvas* ccom = new TCanvas("ccom","ccom",1200,600);
    ccom->Divide(4,2,0,0);
    for(int icent=2;icent<9;icent++){
        if(icent<5) ccom->cd(icent-1);
        else ccom->cd(icent);
        gv1Corrected[0][icent]->SetMaximum(0.18);
        gv1Corrected[0][icent]->SetMinimum(-0.18);
        gv1Corrected[0][icent]->SetMarkerColor(1);
        gv1Corrected[0][icent]->SetLineColor(1);
        gv1Corrected[0][icent]->GetXaxis()->SetTitle("#eta");
        gv1Corrected[0][icent]->GetYaxis()->SetTitle("v_{1}");
        gv1Corrected[0][icent]->Draw("AP");
        gv1Corrected[1][icent]->SetMarkerColor(3);
        gv1Corrected[1][icent]->SetLineColor(3);
        gv1Corrected[1][icent]->Draw("P");
        mv1measured[icent]->Draw("P");
        /*
        for(int ew=0;ew<2;ew++){
            Gurqmdv1[icent][ew]->Draw("P");
        }
        mCosphi[icent]->Draw("same");
        */
        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        l4->Draw("same");
    }
    ccom->cd(8);
    lcom->Draw("same");
    ccom->cd(7);
    leurqmd->Draw("same");
    ccom->SaveAs(OutputPDFName);
    
    TCanvas* curqmd = new TCanvas("curqmd","curqmd",1200,600);
    curqmd->Divide(4,2,0,0);
    for(int icent=2;icent<9;icent++){
        if(icent<5) curqmd->cd(icent-1);
        else curqmd->cd(icent);
        gv1Corrected[0][icent]->SetMaximum(0.18);
        gv1Corrected[0][icent]->SetMinimum(-0.18);
        gv1Corrected[0][icent]->GetXaxis()->SetTitle("#eta");
        gv1Corrected[0][icent]->GetYaxis()->SetTitle("v_{1}");
        gv1Corrected[0][icent]->Draw("AP");
        gv1Corrected[1][icent]->Draw("P");
        for(int ew=0;ew<2;ew++){
            Gurqmdv1[icent][ew]->Draw("P");
        }
        mCosphi[icent]->Draw("same");
        l1->Draw("same");
        l2->Draw("same");
        l3->Draw("same");
        l4->Draw("same");
    }
    curqmd->cd(8);
    leurqmd->Draw("same");
    curqmd->SaveAs(OutputPDFName);

    //combine centrality to compare with PHOBOS 
    double x040[2][16];
    double y040[2][16];
    double xerr040[2][16];
    double yerr040[2][16];

    for(int ich=0;ich<2;ich++){
        for(int iring=0;iring<16;iring++){

            double tmp_sum=0.0;
            double tmp_sum_err=0.0;
            
            for(int icent=8;icent>3;icent--){
                if(icent>6){
                    tmp_sum+=(-yv1InNew[ich][icent][iring]+yv1InNew[ich][icent][31-iring])/16.0;
                    tmp_sum_err+=pow(yerrv1InNew[ich][icent][iring]/16.0,2)+pow(yerrv1InNew[ich][icent][31-iring]/16.0,2);
                }
                else{
                    tmp_sum+=(-yv1InNew[ich][icent][iring]+yv1InNew[ich][icent][31-iring])/8.0;
                    tmp_sum_err+=pow(yerrv1InNew[ich][icent][iring]/8.0,2)+pow(yerrv1InNew[ich][icent][31-iring]/8.0,2);

                }
            }
            y040[ich][iring]=tmp_sum;
            yerr040[ich][iring]=sqrt(tmp_sum_err);
            xerr040[ich][iring]=0.0;
            x040[ich][iring]=xv1InNew[0][8][31-iring]-3.4;
        }

    }
    
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
    
    //v1(TPC) from thr FlowHists_v1mike_run7.root
    double xTPCcom_040[10]={-4.300000,-4.100000,-3.900000,-3.700000,-3.500000,-3.300000,-3.100000,-2.900000,-2.700000,-2.500000};
    double yTPCcom_040[10]={0.008201,0.006361,0.004396,0.002575,0.000821,-0.000989,-0.002790,-0.004515,-0.006401,-0.008346};
    double xerrTPCcom_040[10]={0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000};
    double yerrTPCcom_040[10]={0.000028,0.000027,0.000027,0.000027,0.000027,0.000027,0.000027,0.000027,0.000027,0.000027};

    TCanvas* c1 = new TCanvas("c1","c1",800,600);

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
    TGraphErrors* gr6;
    TGraphErrors* gr7;

    tbg = new TH1D("t1","t1",35, -6.0,3.0);
    tbg->GetYaxis()->SetRangeUser(-0.1,0.2);
    tbg->SetStats(kFALSE);
    tbg->SetTitle("");
    tbg->GetYaxis()->SetTitle("v_{1}");
    tbg->GetXaxis()->SetTitle("#eta-y_{beam}");
    tbg->GetYaxis()->SetLabelSize(0.06);
    tbg->GetXaxis()->SetLabelSize(0.06);
    tbg->GetYaxis()->SetTitleSize(0.06);
    tbg->GetXaxis()->SetTitleSize(0.06);
    tbg->GetYaxis()->SetTitleOffset(0.8);
    tbg->GetXaxis()->SetTitleOffset(0.8);
    tbg->Draw("");

    gr1 = new TGraphErrors(9,x1,y1,xerr1,yerr1);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(1.5);
    gr1->SetMarkerColor(1);
    gr1->SetLineColor(1);
    gr1->Draw("same,P");
    
    gr2 = new TGraphErrors(9,x2,y2,xerr2,yerr2);
    gr2->SetMarkerStyle(23);
    gr2->SetMarkerSize(1.5);
    gr2->SetMarkerColor(8);
    gr2->SetLineColor(8);
    gr2->Draw("same,P");
    cout<<"h1"<<endl; 

    gr3 = new TGraphErrors(9,x3,y3,xerr3,yerr3);
    gr3->SetMarkerStyle(21);
    gr3->SetMarkerSize(1.5);
    gr3->SetMarkerColor(4);
    gr3->SetLineColor(4);
    gr3->Draw("same,P");
    
    gr4 = new TGraphErrors(9,x4,y4,xerr4,yerr4);
    gr4->SetMarkerStyle(26);
    gr4->SetMarkerSize(1.5);
    gr4->SetMarkerColor(2);
    gr4->SetLineColor(2);
    gr4->Draw("same,P");
    
    gr5 = new TGraphErrors(10,xTPCcom_040,yTPCcom_040,xerrTPCcom_040,yerrTPCcom_040);
    gr5->SetMarkerStyle(29);
    gr5->SetMarkerSize(1.5);
    gr5->SetMarkerColor(6);
    gr5->SetLineColor(6);
    gr5->Draw("same,P");

    gr6 = new TGraphErrors(16,x040[0],y040[0],xerr040[0],yerr040[0]);
    gr6->SetMarkerStyle(29);
    gr6->SetMarkerSize(1.5);
    gr6->SetMarkerColor(6);
    gr6->SetLineColor(6);
    gr6->Draw("same,P");
    
    gr7 = new TGraphErrors(16,x040[1],y040[1],xerr040[1],yerr040[1]);
    gr7->SetMarkerStyle(29);
    gr7->SetMarkerSize(1.5);
    gr7->SetMarkerColor(7);
    gr7->SetLineColor(7);
    gr7->Draw("same,P");
    
    TLegend* le1 = new TLegend(0.21,0.58,0.4,0.89);
    le1->AddEntry(gr1,"PHOBOS 19GeV","P");
    le1->AddEntry(gr2,"PHOBOS 62GeV","P");
    le1->AddEntry(gr3,"PHOBOS 130GeV","P");
    le1->AddEntry(gr4,"PHOBOS 200GeV","P");
    le1->AddEntry(gr6,"STAR 27GeV all","P");
    le1->AddEntry(gr7,"STAR 27GeV charged","P");
    le1->SetBorderSize(0);
    le1->SetTextSize(0.03);
    le1->Draw("same");
    c1->SaveAs(OutputPDFName);

    TitlePage->SaveAs(Form("%s%s",OutputPDFName.Data(),"]"));
    
    TString OutputROOTName = "/Users/xiaoyuliu/Library/Mobile Documents/com~apple~CloudDocs/v1/v1_Mike/Output/run26/GEANT/root/"+mOutputFileName+".root";
    TFile *foutput = new TFile(OutputROOTName.Data(),"RECREATE");
    for(int ich=0;ich<2;ich++){
        for(int icent=2;icent<9;icent++){
            gv1Corrected[ich][icent]->SetName(Form("v1CorrectedCh%dCent%d",ich,icent));
            gv1Corrected[ich][icent]->Write(); 
        }
    }
    gr6->SetName("v1_vs_eta_minus_ybeam_all_040");
    gr7->SetName("v1_vs_eta_minus_ybeam_charged_040");
    gr6->Write();
    gr7->Write();
    foutput->Write();
    foutput->Close();
}