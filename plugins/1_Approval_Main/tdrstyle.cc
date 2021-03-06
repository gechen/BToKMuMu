#include "TStyle.h"
#include "TPad.h"

// tdrGrid: Turns the grid lines on (true) or off (false)

// void tdrGrid(bool gridOn) {
//   tdrStyle->SetPadGridX(gridOn);
//   tdrStyle->SetPadGridY(gridOn);
// }

// fixOverlay: Redraws the axis

void fixOverlay() {
  gPad->RedrawAxis();
}

void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(800); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

////// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(3);   //// 3   /////  3 for q2 bins 2D projection plots

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
//  tdrStyle->SetErrorMarker(20);
//  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  tdrStyle->SetMarkerSize(1.3);

//For the fit/function:
  tdrStyle->SetOptFit(0);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

////// For the statistics box:
////  tdrStyle->SetOptFile(0);
////  tdrStyle->SetOptStat(1); // To display the mean and RMS:   SetOptStat("mr");
////  tdrStyle->SetStatColor(kWhite);
////  tdrStyle->SetStatFont(42);
////  tdrStyle->SetStatFontSize(0.08);
////  tdrStyle->SetStatTextColor(1);
////  tdrStyle->SetStatFormat("6.4g");
////  tdrStyle->SetStatBorderSize(2);
////  tdrStyle->SetStatH(0.2);
////  tdrStyle->SetStatW(0.3);
////  // tdrStyle->SetStatStyle(Style_t style = 1001);
////  //tdrStyle->SetStatX(Float_t x = 0.5);
////  //tdrStyle->SetStatY(Float_t y = 0.5);
////
// Margins:
  tdrStyle->SetPadTopMargin(0.09);    //// 0.09   /////   0.09 for q2 bins 2D projection plots
  tdrStyle->SetPadBottomMargin(0.14); //// 0.14   /////   0.19 for q2 bins 2D projection plots
  tdrStyle->SetPadLeftMargin(0.14);   //// 0.14   /////   0.17 for q2 bins 2D projection plots
  tdrStyle->SetPadRightMargin(0.03);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(22);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(22, "XYZ");  // 22
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.13);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(22, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.06, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.04, "XYZ");  //// 0.04  //// 0.04 for for q2 bins 2D projection plots
  tdrStyle->SetNdivisions(505, "XYZ");   //// 505 for efficiency  /////  505  for for q2 bins 2D projection plots
  tdrStyle->SetPadTickX(0);  // 1 -->  To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(0);
  tdrStyle->SetLineWidth(3);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}
