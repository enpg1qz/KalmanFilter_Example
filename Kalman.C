#define max(a,b) (a>b?a:b)

void ComputeRange(Double_t &xmin, Double_t &xmax,Double_t* fX,Int_t fNpoints) const
{
    if(fNpoints <= 0) {
       xmin = xmax = 0;
       return;
    }
    xmin = xmax = fX[0];
    for(Int_t i = 1; i < fNpoints; i++) {
       if(fX[i] < xmin) xmin = fX[i];
       if(fX[i] > xmax) xmax = fX[i];
    }
    //cout<<"min: "<<xmin<<"  max: "<<xmax<<endl;
}


/////https://root.cern.ch/doc/master/TMath_8h.html#a9dc77c64d7b6f9cd8b78e2fd74d63a16
/////https://root.cern.ch/doc/master/namespaceTMath.html#af3ad3bf04aec8eba121523422acdd46d
/////https://root.cern.ch/doc/master/classTComplex.html#a95db81d98b1655255d2f34be856028c7
void FindNiceNdivisions(Int_t &n,Double_t &xmin,Double_t &xmax,Double_t &ymin,Double_t &ymax){
    Double_t L = TMath::Max(xmax-xmin,ymax-ymin)*1.35;
    Int_t m = ceil(TMath::Log10(L))-2;
    Double_t x = L/pow(10,m+1);  //L=x*10^(m+1)  1<x<=10
    if(x>5){
        n = ceil(x);
        cout<<"m: "<<m<<"  n:  "<<n<<"  --1"<<endl;
        //cout<<n<<"  1"<<endl;
        L = 5*n*2*pow(10,m);   /////n1*10^(m+1)~n2*10^(m+1)
        Int_t nx0 = (floor(xmax/pow(10,m+1))+1)-(ceil(xmin/pow(10,m+1))-1);
        if(nx0>n){cout<<"FindNiceNdivisions error: nx0>n"<<endl;}
        Int_t nx_1 = (ceil(xmin/pow(10,m+1))-1)-floor((n-nx0)/2);
        if(((nx_1+n)*pow(10,m+1)-xmax-1)>(xmin-nx_1*pow(10,m+1))){
            nx_1--;
        }
        xmin = nx_1*pow(10,m+1);
        //Double_t xmin_temp = floor((xmin+xmax-L)*0.5/pow(10,m+1))*pow(10,m+1);
        ////if((xmin>xmin_temp+1)&&((xmax>xmin_temp+L)||((xmin-xmin_temp)>(xmin_temp+L-xmax)))){
        //if((xmin>xmin_temp+1)&&(xmax>xmin_temp+L)){
        //    xmin = xmin_temp+1;
        //}
        //else{xmin=xmin_temp;}
        //xmin = floor(xmin/pow(10,m+1))*pow(10,m+1);
        //xmin = (ceil(xmin/pow(10,m+1))-1)*pow(10,m+1);
        xmax = xmin+L;
        //Double_t ymin_temp = floor((ymin+ymax-L)*0.5/pow(10,m+1))*pow(10,m+1);
        Int_t ny0 = (floor(ymax/pow(10,m+1))+1)-(ceil(ymin/pow(10,m+1))-1);
        if(ny0>n){cout<<"FindNiceNdivisions error: ny0>n"<<endl;}
        Int_t ny_1 = ((ceil(ymin/pow(10,m+1))-1)-floor((n-ny0)/2));
        if(((ny_1+n)*pow(10,m+1)-ymax-1)>(ymin-ny_1*pow(10,m+1))){
            ny_1--;
        }
        ymin = ny_1*pow(10,m+1);
        ymax = ymin+L;
        n = 500+n;
    }
    else{
        if(x<3){x = 3;}
        //if(x<2.5){x = 2.5;}
        //if(x<1.5){x = 1.5;}
        n = ceil(x*2);
        cout<<"m: "<<m<<"  n:  "<<n<<"  --1"<<endl;
        //cout<<n<<"  2"<<endl;
        L = 5*n*pow(10,m);   /////n1*5*10^m~n2*5*10^m
        Int_t nx0 = (floor(xmax/(5*pow(10,m)))+1)-(ceil(xmin/(5*pow(10,m)))-1);
        if(nx0>n){cout<<"FindNiceNdivisions error: nx0>n"<<endl;}
        Int_t nx_1 = ((ceil(xmin/(5*pow(10,m)))-1)-floor((n-nx0)/2));
        if(((nx_1+n)*(5*pow(10,m))-xmax-1)>(xmin-nx_1*(5*pow(10,m)))){
            nx_1--;
        }
        xmin = nx_1*(5*pow(10,m));
        //Double_t xmin_temp=floor((xmin+xmax-L)*0.5/(5*pow(10,m)))*(5*pow(10,m));
        ////if((xmin>xmin_temp+1)&&((xmax>xmin_temp+L)||((xmin-xmin_temp)>(xmin_temp+L-xmax)))){
        //if((xmin>xmin_temp+1)&&(xmax>xmin_temp+L)){
        //    xmin = xmin_temp+1;
        //}
        //else{xmin=xmin_temp;}
        //xmin=(ceil(xmin/(5*pow(10,m)))-1)*(5*pow(10,m));
        xmax = xmin+L;
        //Double_t ymin_temp=floor((ymin+ymax-L)*0.5/(5*pow(10,m)))*(5*pow(10,m));
        //if((ymin>ymin_temp+1)&&(ymax>ymin_temp+L)){
        //    ymin = ymin_temp+1;
        //}
        //else{ymin=ymin_temp;}
        //ymin=(ceil(ymin/(5*pow(10,m)))-1)*(5*pow(10,m));
        Int_t ny0 = (floor(ymax/(5*pow(10,m)))+1)-(ceil(ymin/(5*pow(10,m)))-1);
        if(ny0>n){cout<<"FindNiceNdivisions error: ny0>n"<<endl;}
        Int_t ny_1 = ((ceil(ymin/(5*pow(10,m)))-1)-floor((n-ny0)/2));
        if(((ny_1+n)*(5*pow(10,m))-ymax-1)>(ymin-ny_1*(5*pow(10,m)))){
            ny_1--;
        }
        ymin = ny_1*(5*pow(10,m));
        ymax = ymin+L;
        n = 500+n;
    }
}

void SetTMatrixDRow(TMatrixD& A,Int_t row_index,Double_t* A_row,Int_t A_rowUpb){
    if(A.GetNcols()==A_rowUpb){
        TVectorD A_row_temp(A.GetNcols(),A_row);
        TMatrixDRow(A,row_index) = A_row_temp;
    }
    else{
        cout<<A.GetNcols()<<"  !=  "<<A_rowUpb<<endl;
        cout<<"error: A.GetNcols()!=A_rowUpb"<<endl;
        //cout<<"A.GetNcols()!=sizeof(A_row)/sizeof(Double_t)"<<endl;
    }
    return;
}

void PrintTMatrixD(TMatrixD& A,ofstream& outFile){
    for(int m=0;m<A.GetNrows();m++){
        for(int n=0;n<A.GetNcols();n++){
            outFile<<A(m,n)<<"    ";
        }
        outFile<<endl;
    }
}

void PrintTVectorD(TVectorD& A,ofstream& outFile){
    for(int m=0;m<A.GetNrows();m++){
        outFile<<A(m)<<"    ";
    }
    outFile<<endl;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

Double_t omega = 1;
Double_t omega_sigma=0.1;
Double_t eta_sigma = 0.05;
Double_t delta_sigma = 0.1;
Double_t EX_0[4]={0,0,1,2};
Double_t X_0_sigma = 0.1;
Double_t T=0.1;
Int_t N=100;

Double_t t(Int_t k){
    return T*k;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

TMatrixD Phi_Real(Int_t k){  //实际粒子运动使用的Phi ，给omega加上一个随机误差（来模拟磁场的不均匀性）
    Double_t omega_Real = omega+gRandom->Gaus(0,omega_sigma);
    Double_t Phi_0[4] = {1,0,sin(omega_Real*(t(k+1)-t(k)))/omega_Real,(1-cos(omega_Real*(t(k+1)-t(k))))/omega_Real};
    Double_t Phi_1[4] = {0,1,-(1-cos(omega_Real*(t(k+1)-t(k))))/omega_Real,sin(omega_Real*(t(k+1)-t(k)))/omega_Real};
    Double_t Phi_2[4] = {0,0,cos(omega_Real*(t(k+1)-t(k))),sin(omega_Real*(t(k+1)-t(k)))};
    Double_t Phi_3[4] = {0,0,-sin(omega_Real*(t(k+1)-t(k))),cos(omega_Real*(t(k+1)-t(k)))};
    TMatrixD Phi_44(4,4);
    SetTMatrixDRow(Phi_44,0,Phi_0,4);
    SetTMatrixDRow(Phi_44,1,Phi_1,4);
    SetTMatrixDRow(Phi_44,2,Phi_2,4);
    SetTMatrixDRow(Phi_44,3,Phi_3,4);
    return Phi_44;
}

TMatrixD Phi_Filter(Int_t k){  //卡尔曼滤波时估计的Phi
    Double_t Phi_0[4] = {1,0,sin(omega*(t(k+1)-t(k)))/omega,(1-cos(omega*(t(k+1)-t(k))))/omega};
    Double_t Phi_1[4] = {0,1,-(1-cos(omega*(t(k+1)-t(k))))/omega,sin(omega*(t(k+1)-t(k)))/omega};
    Double_t Phi_2[4] = {0,0,cos(omega*(t(k+1)-t(k))),sin(omega*(t(k+1)-t(k)))};
    Double_t Phi_3[4] = {0,0,-sin(omega*(t(k+1)-t(k))),cos(omega*(t(k+1)-t(k)))};
    TMatrixD Phi_44(4,4);
    SetTMatrixDRow(Phi_44,0,Phi_0,4);
    SetTMatrixDRow(Phi_44,1,Phi_1,4);
    SetTMatrixDRow(Phi_44,2,Phi_2,4);
    SetTMatrixDRow(Phi_44,3,Phi_3,4);
    return Phi_44;
}



///Double_t TRandom::Gaus(Double_t mean = 0,Double_t sigma = 1)

TVectorD eta(Int_t k){
    TVectorD eta_k(4);
    Double_t eta_0[4]={0,0,0,0};
    eta_0[2] = gRandom->Gaus(0,eta_sigma);
    eta_0[3] = gRandom->Gaus(0,eta_sigma);
    eta_k.SetElements(eta_0);
    return eta_k;
}

TVectorD get_X0(){
    TVectorD X0(4);
    X0.SetElements(EX_0);
    X0[2] = X0[2]+gRandom->Gaus(0,X_0_sigma);
    X0[3] = X0[3]+gRandom->Gaus(0,X_0_sigma);
    return X0;
}

TMatrixD get_P00(){
    TMatrixD P00(4,4);
    P00.UnitMatrix();
    P00(0,0)=0;
    P00(1,1)=0;
    P00(2,2)=X_0_sigma*X_0_sigma;
    P00(3,3)=X_0_sigma*X_0_sigma;
    return P00;
}


TMatrixD K(Int_t k,TMatrixD **P){
    TMatrixD Kk(4,2);
    if(k<=0){
        Kk = Phi_Real(0)*P[0][0]*H(0).T()*(H(0)*P[0][0]*H(0).T()+R(0)).Invert();
    }
    else{
        Kk = Phi_Real(k)*P[k][k-1]*H(k).T()*(H(k)*P[k][k-1]*H(k).T()+R(k)).Invert();
    }
    return Kk;
}

TMatrixD Q(Int_t k){
    TMatrixD Q_k(4,4);
    Q_k.UnitMatrix();
    Q_k(0,0)=0;
    Q_k(1,1)=0;
    Q_k(2,2)=eta_sigma*eta_sigma;
    Q_k(3,3)=eta_sigma*eta_sigma;
    return Q_k;
}
TMatrixD R(Int_t k){
    TMatrixD R_k(2,2);
    R_k.UnitMatrix();
    R_k(0,0)=delta_sigma*delta_sigma;
    R_k(1,1)=delta_sigma*delta_sigma;
    return R_k;
}

TMatrixD H(Int_t k){
    Double_t H_0[4] = {1,0,0,0};
    Double_t H_1[4] = {0,1,0,0};
    TMatrixD H_42(2,4);
    SetTMatrixDRow(H_42,0,H_0,4);
    SetTMatrixDRow(H_42,1,H_1,4);
    return H_42;
}

TVectorD delta(Int_t k){
    TVectorD delta_k(2);
    Double_t delta_0[2]={0,0};
    delta_0[0] = gRandom->Gaus(0,delta_sigma);
    delta_0[1] = gRandom->Gaus(0,delta_sigma);
    delta_k.SetElements(delta_0);
    return delta_k;
}

/*
TVectorD X(Int_t k){
    TVectorD X_k(4);
    if(k==0){
        Double_t X_0[4]={0,0,1,1};
        X_k.SetElements(X_0);
    }
    else if(k>0){
        X_k = Phi_Real(k-1)*X(k-1)+w(k-1);
    }
    else{
        cout<<"error: k<0"<<endl;
    }
    return X_k;
}
*/


//Int_t N=20;

void Kalman(){
    ofstream outFile;
    //outFile.open("kalman_log.C",ios::app|ios::out);
    outFile.open("kalman_log.C",ios::out);

    delete gRandom;
    gRandom = new TRandom3(0); //seed=0
    ///https://root.cern.ch/doc/master/classTRandom.html
    ///https://root.cern.ch/doc/master/classTUUID.html
    //Int_t N=50;
    TVectorD* X = new TVectorD[N];
    X[0].ResizeTo(4);
    X[0] = get_X0();
    for(Int_t k=1;k<N;k++){
        X[k].ResizeTo(4);
        X[k]=Phi_Real(k-1)*X[k-1]+eta(k-1);
    }
    Double_t* X_x = new Double_t[N];
    Double_t* X_y = new Double_t[N];
    outFile <<"//%%--X[k]--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//"<<endl;
    for(Int_t k=0;k<N;k++){
        X_x[k] = X[k][0];
        X_y[k] = X[k][1];
        //cout<<k<<"   "<<x[k]/y[k]<<endl;
        outFile <<X[k][0]<<"    "<<X[k][1]<<"    "<<X[k][2]<<"    "<<X[k][3]<<endl;
    }

    TVectorD* Z = new TVectorD[N];
    for(Int_t k=0;k<N;k++){
        Z[k].ResizeTo(2);
        Z[k]=H(k)*X[k]+delta(k);
    }
    
    outFile <<"//%%--Z[k]--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//"<<endl;

    Double_t* Z_x = new Double_t[N];
    Double_t* Z_y = new Double_t[N];
    for(Int_t k=0;k<N;k++){
        Z_x[k] = Z[k][0];
        Z_y[k] = Z[k][1];
        //cout<<k<<"   "<<x[k]/y[k]<<endl;
        outFile<<Z[k][0]<<"    "<<Z[k][1]<<endl;
    }

    auto c1 = new TCanvas("c1","c1",0,0,620,650);
    c1->SetCanvasSize(600, 600);
    c1->SetFixedAspectRatio(kTRUE);
    auto grX = new TGraph(N,X_x,X_y);
    //grX->SetLineColor(2);
    //grX->SetLineWidth(4);
    grX->SetMarkerColor(4);
    grX->SetMarkerStyle(6);
    grX->SetMarkerSize(1);
    grX->SetTitle("Kalman Filter");
    grX->GetXaxis()->SetTitle("x");
    grX->GetYaxis()->SetTitle("y");
    //grX->Draw("ACP");
    grX->Draw("ALP");
    //auto xymin = TMath::Min(gPad->GetXmin(),gPad->GetUymax());
    /*
    Double_t xymin = TMath::Min(grX->GetXaxis()->GetXmin(),grX->GetYaxis()->GetXmin());
    Double_t xymax = TMath::Max(grX->GetXaxis()->GetXmax(),grX->GetYaxis()->GetXmax());
    */
    Double_t xmin,xmax,ymin,ymax,xysize;
    ComputeRange(xmin,xmax,X_x,N);
    ComputeRange(ymin,ymax,X_y,N);
    outFile <<"xmin"<<"    "<<"xmax"<<"    "<<"ymin"<<"    "<<"ymax"<<endl;
    outFile <<xmin<<"    "<<xmax<<"    "<<ymin<<"    "<<ymax<<endl;
    Int_t Ndiv;
    FindNiceNdivisions(Ndiv,xmin,xmax,ymin,ymax);
    grX->GetXaxis()->SetNdivisions(Ndiv,kFALSE);
    grX->GetYaxis()->SetNdivisions(Ndiv,kFALSE);
    grX->GetXaxis()->SetLimits(xmin,xmax);
    grX->SetMinimum(ymin);
    grX->SetMaximum(ymax);
    //grX->GetHistogram()->SetMinimum(-20.);
    //grX->GetHistogram()->SetMaximum(20.);
    c1->Modified();

    auto grZ = new TGraph(N,Z_x,Z_y);
    grZ->SetMarkerColor(6);
    grZ->SetMarkerStyle(7);
    grZ->SetMarkerSize(1);
    //grZ->SetTitle("Kalman Filter");
    //grZ->GetXaxis()->SetTitle("x");
    //grZ->GetYaxis()->SetTitle("y");
    //grZ->Draw("CP");
    grZ->Draw("LP");
    c1->Modified();

    // Do Kalman Filter
    TMatrixD **P = new TMatrixD*[N]; //Row
    for(Int_t m=0;m<N;m++){
        P[m] = new TMatrixD[N]; //Col
    }
    P[0][0].ResizeTo(4,4);
    P[0][0]=get_P00();
    P[1][0].ResizeTo(4,4);
    P[1][0]=Phi_Filter(0)*P[0][0]*Phi_Filter(0).T()-Phi_Filter(0)*P[0][0]*H(0).T()*(H(0)*P[0][0]*H(0).T()+R(0)).Invert()*H(0)*P[0][0]*Phi_Filter(0).T()+Q(0);
    for(Int_t k=1;k<N-1;k++){
        P[k+1][k].ResizeTo(4,4);
        P[k+1][k]=Phi_Filter(k)*P[k][k-1]*Phi_Filter(k).T()-Phi_Filter(k)*P[k][k-1]*H(k).T()*(H(k)*P[k][k-1]*H(k).T()+R(k)).Invert()*H(k)*P[k][k-1]*Phi_Filter(k).T()+Q(k);
        cout<<k<<endl;
    }

    TVectorD **X_filter = new TVectorD*[N]; //Row
    for(Int_t m=0;m<N;m++){
        X_filter[m] = new TVectorD[N]; //Col
    }
    X_filter[0][0].ResizeTo(4);
    X_filter[0][0].SetElements(EX_0);
    X_filter[1][0].ResizeTo(4);
    X_filter[1][0]=Phi_Filter(0)*X_filter[0][0]+K(0,P)*(Z[0]-H(0)*X_filter[0][0]);
    for(Int_t k=1;k<N-1;k++){
        X_filter[k+1][k].ResizeTo(4);
        X_filter[k+1][k]=Phi_Filter(k)*X_filter[k][k-1]+K(k,P)*(Z[k]-H(k)*X_filter[k][k-1]);
        cout<<k<<endl;
    }

    outFile<<"//%%--X_filter[k]--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//"<<endl;

    Double_t* X_filter_x = new Double_t[N];
    Double_t* X_filter_y = new Double_t[N];
    X_filter_x[0]=X_filter[0][0].GetMatrixArray()[0];
    X_filter_y[0]=X_filter[0][0].GetMatrixArray()[1];
    for(Int_t k=1;k<N;k++){
        X_filter_x[k] = X_filter[k][k-1].GetMatrixArray()[0];
        X_filter_y[k] = X_filter[k][k-1].GetMatrixArray()[1];
        //cout<<k<<"   "<<x[k]/y[k]<<endl;
        //cout<<X_filter[k][k-1](0)<<endl;
        outFile<<X_filter[k][k-1].GetMatrixArray()[0]<<"    "<<X_filter[k][k-1].GetMatrixArray()[1]<<"    "<<X_filter[k][k-1].GetMatrixArray()[2]<<"    "<<X_filter[k][k-1].GetMatrixArray()[3]<<endl;
    }

    auto grX_F = new TGraph(N,X_filter_x,X_filter_y);
    //TColor * mycolor = TColor();
    //mycolor->SetRGB(0.5,0.5,0.5);// r,g,b  0~1
    grX_F->SetMarkerColor(2);
    grX_F->SetMarkerStyle(7);
    grX_F->SetMarkerSize(1);
    //grX_F->SetTitle("Kalman Filter");
    //grX_F->GetXaxis()->SetTitle("x");
    //grX_F->GetYaxis()->SetTitle("y");
    //grX_F->Draw("CP");
    grX_F->Draw("LP");
    c1->Modified();

    TLegend leg(0.7,0.8,0.9,0.9);
    leg.SetFillColor(0);
    //cout<<"hhh"<<endl;
    leg.AddEntry(grX,"Real Point","p");
    leg.AddEntry(grZ,"Measurement Point","p");
    leg.AddEntry(grX_F,"Estimate Point","p");
    //leg.AddEntry("","LISE","p");
    leg.DrawClone("Same");


    outFile <<"//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//"<<endl;

    outFile<<"P[0][0]:"<<endl;
    PrintTMatrixD(P[0][0],outFile);
    outFile<<"P[1][0]:"<<endl;
    PrintTMatrixD(P[1][0],outFile);
    for(int k =1;k<TMath::Min(N-1,2);k++){
        outFile<<"-----------------------------"<<endl;
        outFile<<"k= "<<k<<endl;
        outFile<<"P[k][k-1]:"<<endl;
        PrintTMatrixD(P[k][k-1],outFile);
        outFile<<"Phi_Filter(k):"<<endl;
        PrintTMatrixD(Phi_Filter(k),outFile);
        outFile<<"R(k):"<<endl;
        PrintTMatrixD(R(k),outFile);
        outFile<<"H(k):"<<endl;
        PrintTMatrixD(H(k),outFile);
        outFile<<"Q(k):"<<endl;
        PrintTMatrixD(Q(k),outFile);
        outFile<<"P[k+1][k]:"<<endl;
        PrintTMatrixD(P[k+1][k],outFile);

        outFile<<"X_filter[k][k-1]:"<<endl;
        PrintTVectorD(X_filter[k][k-1],outFile);
        outFile<<"K(k,P):"<<endl;
        PrintTMatrixD(K(k,P),outFile);
        outFile<<"Z[k]:"<<endl;
        PrintTVectorD(Z[k],outFile);
        outFile<<"X_filter[k+1][k]:"<<endl;
        PrintTVectorD(X_filter[k+1][k],outFile);
    }




    c1->SaveAs("kalman.root");
    outFile.close();


    delete[] X_x;
    delete[] X_y;
    delete[] X;
    delete[] Z_x;
    delete[] Z_y;
    delete[] Z;

    for(Int_t m=0; m<N;m++){
        delete[] P[m];
        delete[] X_filter[m];
    }
    delete[] P;
    delete[] X_filter;
}