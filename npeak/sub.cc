#include "sub.h"

void FillTriLo(double value_bin, double area_bin,
               const Vect2d* const vect_c0,
               const Vect2d* const vect_c1,
               const Vect2d* const vect_c2,
               double y_lo, double y_up,
               HistDataNerr1d* const hd1d_proj)
{
    int ibin_c0 = hd1d_proj->GetHi1d()->GetIbin( vect_c0->GetPosX() );
    int ibin_c1 = hd1d_proj->GetHi1d()->GetIbin( vect_c1->GetPosX() );
    
    if(ibin_c0 == ibin_c1){
        double y_cross0 = vect_c1->GetPosY();
        double y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       vect_c1->GetPosX());
        double area3 = fabs(y_cross0 - y_cross1) *
            (vect_c1->GetPosX() - vect_c0->GetPosX()) / 2.0;
        hd1d_proj->Fill(vect_c0->GetPosX(), value_bin * area3 / area_bin);
    } else if(ibin_c0 + 1 == ibin_c1){
        double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin_c0);
        double y_cross0 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c1->GetPosX(),
                                       vect_c1->GetPosY(),
                                       binup);
        double y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       binup);
        double area3 = fabs(y_cross0 - y_cross1) *
            (binup - vect_c0->GetPosX()) / 2.0;
        hd1d_proj->Fill(vect_c0->GetPosX(), value_bin * area3 / area_bin);
        double y_cross2 = vect_c1->GetPosY();
        double y_cross3 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       vect_c1->GetPosX());
        double area4 = (fabs(y_cross0 - y_cross1) + fabs(y_cross2 - y_cross3)) *
            (vect_c1->GetPosX() - binup) / 2.0;
        hd1d_proj->Fill(vect_c1->GetPosX(), value_bin * area4 / area_bin);
    } else if(ibin_c0 + 1 < ibin_c1){
        double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin_c0);
        double y_cross0 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c1->GetPosX(),
                                       vect_c1->GetPosY(),
                                       binup);
        double y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       binup);
        double area3 = fabs(y_cross0 - y_cross1) *
            (binup - vect_c0->GetPosX()) / 2.0;
        hd1d_proj->Fill(vect_c0->GetPosX(), value_bin * area3 / area_bin);

        for(int ibin = ibin_c0 + 1; ibin < ibin_c1; ibin ++){
            double binlo = hd1d_proj->GetHi1d()->GetBinLo(ibin);
            double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin);
            double bin_center = hd1d_proj->GetHi1d()->GetBinCenter(ibin);
            double y_cross0 = GetCrossYPos(vect_c0->GetPosX(),
                                           vect_c0->GetPosY(),
                                           vect_c1->GetPosX(),
                                           vect_c1->GetPosY(),
                                           binlo);
            double y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                           vect_c0->GetPosY(),
                                           vect_c2->GetPosX(),
                                           vect_c2->GetPosY(),
                                           binlo);
            double y_cross2 = GetCrossYPos(vect_c0->GetPosX(),
                                           vect_c0->GetPosY(),
                                           vect_c1->GetPosX(),
                                           vect_c1->GetPosY(),
                                           binup);
            double y_cross3 = GetCrossYPos(vect_c0->GetPosX(),
                                           vect_c0->GetPosY(),
                                           vect_c2->GetPosX(),
                                           vect_c2->GetPosY(),
                                           binup);
            double area4 = (fabs(y_cross0 - y_cross1) + fabs(y_cross2 - y_cross3)) *
                (binup - binlo) / 2.0;
            hd1d_proj->Fill(bin_center, value_bin * area4 / area_bin);
        }
        double binlo = hd1d_proj->GetHi1d()->GetBinLo(ibin_c1);
        y_cross0 = GetCrossYPos(vect_c0->GetPosX(),
                                vect_c0->GetPosY(),
                                vect_c1->GetPosX(),
                                vect_c1->GetPosY(),
                                binlo);
        y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                vect_c0->GetPosY(),
                                vect_c2->GetPosX(),
                                vect_c2->GetPosY(),
                                binlo);
        double y_cross2 = vect_c1->GetPosY();
        double y_cross3 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       vect_c1->GetPosX());
        double area4 = (fabs(y_cross0 - y_cross1) + fabs(y_cross2 - y_cross3)) *
            (vect_c1->GetPosX() - binlo) / 2.0;
        hd1d_proj->Fill(vect_c1->GetPosX(), value_bin * area4 / area_bin);
    }
}


void FillTriUp(double value_bin, double area_bin,
               const Vect2d* const vect_c1,
               const Vect2d* const vect_c2,
               const Vect2d* const vect_c3,
               double y_lo, double y_up,
               HistDataNerr1d* const hd1d_proj)
{
    int ibin_c2 = hd1d_proj->GetHi1d()->GetIbin( vect_c2->GetPosX() );
    int ibin_c3 = hd1d_proj->GetHi1d()->GetIbin( vect_c3->GetPosX() );
    
    if(ibin_c2 == ibin_c3){
        double y_cross0 = vect_c2->GetPosY();
        double y_cross1 = GetCrossYPos(vect_c1->GetPosX(),
                                       vect_c1->GetPosY(),
                                       vect_c3->GetPosX(),
                                       vect_c3->GetPosY(),
                                       vect_c2->GetPosX());
        double area3 = fabs(y_cross0 - y_cross1) *
            (vect_c3->GetPosX() - vect_c2->GetPosX()) / 2.0;
        hd1d_proj->Fill(vect_c2->GetPosX(), value_bin * area3 / area_bin);
    } else if(ibin_c2 + 1 == ibin_c3){
        double y_cross0 = vect_c2->GetPosY();
        double y_cross1 = GetCrossYPos(vect_c1->GetPosX(),
                                       vect_c1->GetPosY(),
                                       vect_c3->GetPosX(),
                                       vect_c3->GetPosY(),
                                       vect_c2->GetPosX());
        double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin_c2);
        double y_cross2 = GetCrossYPos(vect_c1->GetPosX(),
                                       vect_c1->GetPosY(),
                                       vect_c3->GetPosX(),
                                       vect_c3->GetPosY(),
                                       binup);
        double y_cross3 = GetCrossYPos(vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       vect_c3->GetPosX(),
                                       vect_c3->GetPosY(),
                                       binup);
        double area4 = (fabs(y_cross0 - y_cross1) + fabs(y_cross2 - y_cross3)) *
            (binup - vect_c2->GetPosX()) / 2.0;
        hd1d_proj->Fill(vect_c2->GetPosX(), value_bin * area4 / area_bin);
        double area3 = fabs(y_cross2 - y_cross3) *
            (vect_c3->GetPosX() - binup) / 2.0;
        hd1d_proj->Fill(vect_c3->GetPosX(), value_bin * area3 / area_bin);
    } else if(ibin_c2 + 1 < ibin_c3){
        double y_cross0 = vect_c2->GetPosY();
        double y_cross1 = GetCrossYPos(vect_c1->GetPosX(),
                                       vect_c1->GetPosY(),
                                       vect_c3->GetPosX(),
                                       vect_c3->GetPosY(),
                                       vect_c2->GetPosX());
        double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin_c2);
        double y_cross2 = GetCrossYPos(vect_c1->GetPosX(),
                                       vect_c1->GetPosY(),
                                       vect_c3->GetPosX(),
                                       vect_c3->GetPosY(),
                                       binup);
        double y_cross3 = GetCrossYPos(vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       vect_c3->GetPosX(),
                                       vect_c3->GetPosY(),
                                       binup);
        double area4 = (fabs(y_cross0 - y_cross1) + fabs(y_cross2 - y_cross3)) *
            (binup - vect_c2->GetPosX()) / 2.0;
        hd1d_proj->Fill(vect_c2->GetPosX(), value_bin * area4 / area_bin);

        for(int ibin = ibin_c2 + 1; ibin < ibin_c3; ibin ++){
            double binlo = hd1d_proj->GetHi1d()->GetBinLo(ibin);
            double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin);
            double bin_center = hd1d_proj->GetHi1d()->GetBinCenter(ibin);
            double y_cross0 = GetCrossYPos(vect_c1->GetPosX(),
                                           vect_c1->GetPosY(),
                                           vect_c3->GetPosX(),
                                           vect_c3->GetPosY(),
                                           binlo);
            double y_cross1 = GetCrossYPos(vect_c2->GetPosX(),
                                           vect_c2->GetPosY(),
                                           vect_c3->GetPosX(),
                                           vect_c3->GetPosY(),
                                           binlo);
            double y_cross2 = GetCrossYPos(vect_c1->GetPosX(),
                                           vect_c1->GetPosY(),
                                           vect_c3->GetPosX(),
                                           vect_c3->GetPosY(),
                                           binup);
            double y_cross3 = GetCrossYPos(vect_c2->GetPosX(),
                                           vect_c2->GetPosY(),
                                           vect_c3->GetPosX(),
                                           vect_c3->GetPosY(),
                                           binup);
            double area4 = (fabs(y_cross0 - y_cross1) + fabs(y_cross2 - y_cross3)) *
                (binup - binlo) / 2.0;
            hd1d_proj->Fill(bin_center, value_bin * area4 / area_bin);
        }
        double binlo = hd1d_proj->GetHi1d()->GetBinLo(ibin_c3);
        y_cross0 = GetCrossYPos(vect_c1->GetPosX(),
                                vect_c1->GetPosY(),
                                vect_c3->GetPosX(),
                                vect_c3->GetPosY(),
                                binlo);
        y_cross1 = GetCrossYPos(vect_c2->GetPosX(),
                                vect_c2->GetPosY(),
                                vect_c3->GetPosX(),
                                vect_c3->GetPosY(),
                                binlo);
        double area3 = fabs(y_cross0 - y_cross1) *
            (vect_c3->GetPosX() - binlo) / 2.0;
        hd1d_proj->Fill(vect_c3->GetPosX(), value_bin * area3 / area_bin);
    }
}

// fill parallelogram
void FillSq(double value_bin, double area_bin,
            const Vect2d* const vect_c0,
            const Vect2d* const vect_c1,
            const Vect2d* const vect_c2,
            const Vect2d* const vect_c3,
            double y_lo, double y_up,
            HistDataNerr1d* const hd1d_proj)
{
    int ibin_c1 = hd1d_proj->GetHi1d()->GetIbin( vect_c1->GetPosX() );
    int ibin_c2 = hd1d_proj->GetHi1d()->GetIbin( vect_c2->GetPosX() );
    
    if(ibin_c1 == ibin_c2){
        double y_cross0 = vect_c1->GetPosY();
        double y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       vect_c1->GetPosX());
        double area3 = fabs(y_cross0 - y_cross1) *
            (vect_c2->GetPosX() - vect_c1->GetPosX());
        hd1d_proj->Fill(vect_c1->GetPosX(), value_bin * area3 / area_bin);
    } else if(ibin_c1 + 1 == ibin_c2){
        double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin_c1);
        double y_cross0 = vect_c1->GetPosY();
        double y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       vect_c1->GetPosX());
        double area4 = fabs(y_cross0 - y_cross1) *
            (binup - vect_c1->GetPosX());
        hd1d_proj->Fill(vect_c1->GetPosX(), value_bin * area4 / area_bin);
        area4 = fabs(y_cross0 - y_cross1) *
            (vect_c2->GetPosX() - binup);
        hd1d_proj->Fill(vect_c2->GetPosX(), value_bin * area4 / area_bin);
    } else if(ibin_c1 + 1 < ibin_c2){
        double y_cross0 = vect_c1->GetPosY();
        double y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                       vect_c0->GetPosY(),
                                       vect_c2->GetPosX(),
                                       vect_c2->GetPosY(),
                                       vect_c1->GetPosX());
        double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin_c1);
        double area4 = fabs(y_cross0 - y_cross1) *
            (binup - vect_c1->GetPosX());
        hd1d_proj->Fill(vect_c1->GetPosX(), value_bin * area4 / area_bin);

        for(int ibin = ibin_c1 + 1; ibin < ibin_c2; ibin ++){
            double binlo = hd1d_proj->GetHi1d()->GetBinLo(ibin);
            double binup = hd1d_proj->GetHi1d()->GetBinUp(ibin);
            double bin_center = hd1d_proj->GetHi1d()->GetBinCenter(ibin);
            double y_cross0 = GetCrossYPos(vect_c1->GetPosX(),
                                           vect_c1->GetPosY(),
                                           vect_c3->GetPosX(),
                                           vect_c3->GetPosY(),
                                           binlo);
            double y_cross1 = GetCrossYPos(vect_c0->GetPosX(),
                                           vect_c0->GetPosY(),
                                           vect_c2->GetPosX(),
                                           vect_c2->GetPosY(),
                                           binlo);
            double area4 = fabs(y_cross0 - y_cross1) * (binup - binlo);
            hd1d_proj->Fill(bin_center, value_bin * area4 / area_bin);
        }
        double binlo = hd1d_proj->GetHi1d()->GetBinLo(ibin_c2);
        y_cross0 = vect_c2->GetPosY();
        y_cross1 = GetCrossYPos(vect_c1->GetPosX(),
                                vect_c1->GetPosY(),
                                vect_c3->GetPosX(),
                                vect_c3->GetPosY(),
                                vect_c2->GetPosX());
        area4 = fabs(y_cross0 - y_cross1) * (vect_c2->GetPosX() - binlo);
        hd1d_proj->Fill(vect_c2->GetPosX(), value_bin * area4 / area_bin);
    }
}



double GetCrossYPos(double xpos1, double ypos1,
                    double xpos2, double ypos2,
                    double xpos)
{
    double ypos = ypos1 + (ypos2 - ypos1) * (xpos - xpos1) / (xpos2 - xpos1);
    return ypos;
}

