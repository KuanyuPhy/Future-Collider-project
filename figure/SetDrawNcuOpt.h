inline void SetRootStyle()
{
    //  For the gstyle statistics:
    //gStyle->SetOptStat(0);
    gStyle->SetOptFit(111111);

}

template <typename T>
void SetDrawNcuOpt(T &h, Color_t Color, string title, string xTitle, string yTitle)
{

    //  For the Global title:
    h->SetTitle(title.c_str());

    //  For the X axis titles & Lable:
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleSize(0.05);
    h->SetXTitle(xTitle.c_str());

    //  For the Y axis titles & Lable:

    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetTitleSize(0.05);
    h->SetYTitle(yTitle.c_str());

    //  For the Line Style
    h->SetLineColor(Color);
    h->SetLineWidth(2);

    // For the statistics box:
    h->SetStats(1);
}
